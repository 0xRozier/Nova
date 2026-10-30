"""Minimal stand-in for the `bittensor` package, for running BoltzGen alone.

The NOVA-vendored BoltzGen imports bittensor solely to log (nine
`bt.logging.{debug,info,error}` calls in boltzgen_wrapper.py and
cli/boltzgen.py). Installing the real package into BoltzGen's pinned
environment risks dependency conflicts for no benefit, so the scoring notebook
drops this file into that interpreter's site-packages instead.

Never use this in a miner or validator: it provides no chain access.
"""
import sys


class _Logging:
    """Prints to stderr, mirroring bittensor.logging's surface."""

    @staticmethod
    def _emit(level: str, *args, **kwargs):
        print(f"[boltzgen:{level}]", *args, file=sys.stderr, flush=True)

    def debug(self, *args, **kwargs):
        self._emit("debug", *args, **kwargs)

    def info(self, *args, **kwargs):
        self._emit("info", *args, **kwargs)

    def warning(self, *args, **kwargs):
        self._emit("warning", *args, **kwargs)

    def error(self, *args, **kwargs):
        self._emit("error", *args, **kwargs)

    def trace(self, *args, **kwargs):
        self._emit("trace", *args, **kwargs)

    def success(self, *args, **kwargs):
        self._emit("success", *args, **kwargs)

    def set_trace(self, *args, **kwargs):
        pass

    def __getattr__(self, name):
        return lambda *args, **kwargs: None


logging = _Logging()
__version__ = "0.0.0+nova-lab-stub"
