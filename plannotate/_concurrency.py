"""Core allocation and bounded concurrency for annotation sources."""

import shlex
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from typing import Any, Callable, TypeVar

Result = TypeVar("Result")
Search = Callable[[str, str, dict[str, Any], bool, int], Result]


def parameters_with_threads(
    parameters: str,
    thread_options: tuple[str, ...],
    thread_option: str,
    threads: int,
) -> str:
    """Replace configured thread flags with the allocated thread count."""
    if threads < 1:
        raise ValueError("threads must be at least 1")

    tokens = iter(shlex.split(parameters))
    filtered = []
    for token in tokens:
        if token in thread_options:
            try:
                next(tokens)
            except StopIteration as exc:
                raise ValueError(f"Thread option {token!r} requires a value") from exc
            continue
        if any(token.startswith(f"{option}=") for option in thread_options):
            continue
        filtered.append(token)
    filtered.extend((thread_option, str(threads)))
    return shlex.join(filtered)


def _lane_allocation(costs: list[float], cores: int, lanes: int) -> list[int]:
    """Give each of ``lanes`` heaviest sources a beefy thread count, others one.

    The ``lanes`` highest-cost sources are run in parallel and split the core
    budget between them (greedily, by projected ``cost / threads`` bottleneck);
    every remaining source is a one-thread filler that runs behind a lane when it
    frees. The heavy allocations sum to ``cores``, so at most ``lanes`` concurrent
    tasks never exceed the budget.
    """
    allocations = [1] * len(costs)
    primaries = sorted(range(len(costs)), key=lambda index: costs[index], reverse=True)
    primaries = primaries[:lanes]
    for _ in range(cores - lanes):
        bottleneck = max(primaries, key=lambda index: costs[index] / allocations[index])
        allocations[bottleneck] += 1
    return allocations


def _estimate_makespan(times: list[float], worker_count: int) -> float:
    """Estimate wall-clock time of list-scheduling tasks (in order) onto workers."""
    worker_free = [0.0] * worker_count
    for duration in times:
        soonest = min(range(worker_count), key=lambda index: worker_free[index])
        worker_free[soonest] += duration
    return max(worker_free)


def plan_concurrency(
    sources: dict[str, dict[str, Any]], cores: int
) -> tuple[int, list[int]]:
    """Choose how many sources run at once and how many threads each one gets.

    Sources run concurrently, so wall-clock time is the slowest lane, not the sum.
    Fewer, beefier lanes let a scalable bottleneck (e.g. an Infernal search) take
    several threads instead of being pinned to one when there are more sources
    than cores; trivial sources just queue behind a lane. We pick the lane count
    that minimizes the projected makespan under near-linear ``cost / threads``
    scaling, never exceeding the core budget. With equal costs and ample cores
    this reduces to one balanced thread per source.
    """
    if cores < 1:
        raise ValueError("cores must be at least 1")
    if not sources:
        return 0, []

    costs = [float(config.get("cost", 1.0)) for config in sources.values()]
    best_makespan = float("inf")
    best_workers = min(cores, len(sources))
    best_allocations = [1] * len(sources)
    # try most lanes first so ties keep maximum concurrency (closest to balanced)
    for lanes in range(min(cores, len(sources)), 0, -1):
        allocations = _lane_allocation(costs, cores, lanes)
        times = [costs[index] / allocations[index] for index in range(len(costs))]
        makespan = _estimate_makespan(times, lanes)
        if makespan < best_makespan - 1e-9:
            best_makespan = makespan
            best_workers = lanes
            best_allocations = allocations
    return best_workers, best_allocations


def run_sources(
    search: Search[Result],
    query: str,
    is_linear: bool,
    sources: dict[str, dict[str, Any]],
    cores: int,
) -> list[Result]:
    """Run annotation sources concurrently and return results in YAML order."""
    if not sources:
        return []

    names = list(sources)
    configs = list(sources.values())
    worker_count, threads = plan_concurrency(sources, cores)
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        return list(
            executor.map(
                search,
                repeat(query),
                names,
                configs,
                repeat(is_linear),
                threads,
            )
        )
