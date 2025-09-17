"""
Performance Optimizer
--------------------
Dynamic resource allocation, memory profiling, and caching for high-throughput sequence processing.
Ensures robust, efficient, and scalable pipeline execution.
"""

import logging
from typing import Any, Callable

try:
    import psutil  # type: ignore
except ImportError:  # pragma: no cover
    psutil = None  # type: ignore

try:
    import joblib  # type: ignore
except ImportError:  # pragma: no cover
    joblib = None  # type: ignore

class PerformanceOptimizerError(Exception):
    pass

class PerformanceOptimizer:
    def __init__(self):
        self.logger = logging.getLogger("PerformanceOptimizer")
        self.memory_limit = (psutil.virtual_memory().total * 0.8) if psutil else 0  # Use up to 80% of available RAM
        self.cache = joblib.Memory(location=".performance_cache", verbose=0) if joblib else None

    def dynamic_resource_allocation(self, dataset_size: int) -> int:
        """
        Dynamically determines optimal number of threads based on dataset size and available resources.
        """
        cpu_count = psutil.cpu_count(logical=False) if psutil else None
        if cpu_count is None:
            cpu_count = 1
        if dataset_size < 1000:
            return max(1, cpu_count // 2)
        elif dataset_size < 10000:
            return max(1, cpu_count - 1)
        else:
            return cpu_count

    def memory_profiler(self, func: Callable, *args, **kwargs) -> Any:
        """
        Profiles memory usage of a function call.
        """
        import tracemalloc
        tracemalloc.start()
        result = func(*args, **kwargs)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.logger.info(f"Memory usage: current={current / 1e6:.2f}MB, peak={peak / 1e6:.2f}MB")
        return result

    def cache_result(self, func: Callable) -> Callable:
        """
        Decorator to cache function results for efficiency.
        """
        if self.cache is None:
            return func
        return self.cache.cache(func)
