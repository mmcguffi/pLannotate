"""Core allocation and bounded concurrency for annotation sources."""

import shlex
from concurrent.futures import ThreadPoolExecutor
from itertools import cycle, islice, repeat
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


def allocate_threads(sources: dict[str, dict[str, Any]], cores: int) -> list[int]:
    """Allocate the core budget across configured annotation sources."""
    if cores < 1:
        raise ValueError("cores must be at least 1")
    if not sources:
        return []

    allocations = [1] * len(sources)
    spare_cores = max(0, cores - len(sources))
    source_indexes = range(len(sources))
    for index in islice(cycle(source_indexes), spare_cores):
        allocations[index] += 1
    return allocations


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
    threads = allocate_threads(sources, cores)
    worker_count = min(cores, len(sources))
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
