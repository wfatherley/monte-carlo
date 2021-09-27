"""gillespy"""
from logging import (
    getLogger,
    handlers,
    StreamHandler,
    ERROR,
    INFO
)
from sys import platform

from .method import Base, Direct, FirstFamily, FirstReaction
from .model import Model
from .util import GillespyException


stream_handler = StreamHandler()
stream_handler.setLevel(ERROR)

if platform == "cygwin":
    system_log_handler = handlers.NTEventLogHandler()
else:
    system_log_handler = handlers.SysLogHandler()
system_log_handler.setLevel(INFO)

logger = getLogger("gillespy")
logger.addHandler(stream_handler)
logger.addHandler(system_log_handler)