import os
import signal

import pytest


DEFAULT_TEST_TIMEOUT_SECONDS = 120
DEFAULT_INTEGRATION_TIMEOUT_SECONDS = 900
TRUE_VALUES = {"1", "true", "yes"}


def _nonnegative_int(config, option_name):
    value = config.getoption(option_name)
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise pytest.UsageError(f"{option_name} must be an integer") from exc
    if parsed < 0:
        raise pytest.UsageError(f"{option_name} must be nonnegative")
    return parsed


def pytest_addoption(parser):
    group = parser.getgroup("plannotate")
    group.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="run tests that require external tools and downloaded databases",
    )
    group.addoption(
        "--strict-annotation-controls",
        action="store_true",
        default=os.environ.get("PLANNOTATE_STRICT_ANNOTATION_CONTROLS", "").lower()
        in TRUE_VALUES,
        help="fail instead of xfail when annotation output differs from controls",
    )
    group.addoption(
        "--test-timeout",
        default=os.environ.get(
            "PLANNOTATE_TEST_TIMEOUT", str(DEFAULT_TEST_TIMEOUT_SECONDS)
        ),
        help="seconds before a unit test fails; use 0 to disable",
    )
    group.addoption(
        "--integration-timeout",
        default=os.environ.get(
            "PLANNOTATE_INTEGRATION_TEST_TIMEOUT",
            str(DEFAULT_INTEGRATION_TIMEOUT_SECONDS),
        ),
        help="seconds before an integration test fails; use 0 to disable",
    )


def pytest_configure(config):
    config.plannotate_test_timeout = _nonnegative_int(config, "--test-timeout")
    config.plannotate_integration_timeout = _nonnegative_int(
        config, "--integration-timeout"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-integration") or config.option.collectonly:
        return
    selected = []
    deselected = []
    for item in items:
        (deselected if "integration" in item.keywords else selected).append(item)
    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = selected


def _timeout_for(item):
    marker = item.get_closest_marker("timeout")
    if marker is not None:
        try:
            return int(marker.args[0])
        except (IndexError, TypeError, ValueError) as exc:
            raise pytest.UsageError(
                f"{item.nodeid} has invalid timeout marker; expected timeout(seconds)"
            ) from exc
    if "integration" in item.keywords:
        return item.config.plannotate_integration_timeout
    return item.config.plannotate_test_timeout


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    timeout = _timeout_for(item)
    if timeout == 0 or not hasattr(signal, "SIGALRM"):
        yield
        return

    def handle_timeout(signum, frame):
        raise TimeoutError(f"{item.nodeid} exceeded {timeout}s test timeout")

    try:
        previous_handler = signal.signal(signal.SIGALRM, handle_timeout)
    except ValueError:
        yield
        return

    signal.alarm(timeout)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)
