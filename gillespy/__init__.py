"""gillespy: stochastic simulation algorithms"""
from logging import getLogger, handlers, StreamHandler, ERROR, INFO
from sys import platform

from .method import Base, Direct, FirstReaction
from .model import Model


stream_handler = StreamHandler()
stream_handler.setLevel(ERROR)

if platform == "cygwin":
    system_log_handler = handlers.NTEventLogHandler()
else:
    system_log_handler = handlers.SysLogHandler()
system_log_handler.setLevel(INFO)

logger = getLogger(__name__)
logger.addHandler(stream_handler)
logger.addHandler(system_log_handler)

__all__ = ["Base", "Direct", "FirstReaction", "Model"]
