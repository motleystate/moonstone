from functools import lru_cache
from logging import getLogger

# Keep track of 10 different messages and then warn again
@lru_cache(10)
def warn_once(logger: getLogger, msg: str):
    """"
    To not repeat the same warning messages more than every n=10 differerent messages.
    NB: only track messages passing through this function. 
    
    by Kound
    src: https://stackoverflow.com/questions/31953272/logging-print-message-only-once
    """
    logger.warning(msg)

def reset_warnings():
    warn_once.cache_clear()

def reset_warnings_decorator(func):
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        reset_warnings()  # Clear warning cache after function execution
        return result
    return wrapper