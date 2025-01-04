from collections import defaultdict, deque
import threading
import time
import requests

class APITracker:
    def __init__(self):
        self.session = requests.Session()
        self.ncbi_call_count = defaultdict(int)
        self.total_calls = 0
        self.all_urls = []
        
    def request(self, method, url, *args, **kwargs):
        # Increment total calls
        self.total_calls += 1
        
        # Store all URLs
        self.all_urls.append(url)
        
        # Track NCBI specific calls
        if "ncbi" in url.lower():
            self.ncbi_call_count[f"{method} {url}"] += 1
            
        return self.session.request(method, url, *args, **kwargs)
    
    def get_stats(self):
        return {
            "total_api_calls": self.total_calls,
            "ncbi_calls": dict(self.ncbi_call_count),
            "ncbi_total_calls": sum(self.ncbi_call_count.values()),
            "all_urls": self.all_urls
        }

tracker = APITracker()

# Simple in-memory rate limiter for NCBI requests: no more than 9 rps
# We'll track timestamps of requests in a deque
ncbi_requests_log = deque()
ncbi_lock = threading.Lock()
# Global counters
total_api_time = 0.0
api_call_count = 0
# Global session (persistent)
session = requests.Session()
session.headers.update({
    "User-Agent": "ligify-agent/1.0",
    # Include any other headers you may need globally
})

# Desired rate limit: 9 requests/second means ~0.111s between requests.
# We'll pick 0.12s to be safe. Adjust if needed.
RATE_LIMIT_INTERVAL = 0.12

def make_request(method, url, **kwargs):
    global total_api_time, api_call_count

    start_time = time.time()

    # Perform the request using the persistent session
    resp = session.request(method, url, **kwargs)

    # If you need retry logic, add it here. For example:
    # retry_count = 0
    # while resp.status_code == 429 and retry_count < 5:
    #     time.sleep(1)
    #     resp = session.request(method, url, **kwargs)
    #     retry_count += 1

    end_time = time.time()
    duration = end_time - start_time
    total_api_time += duration
    api_call_count += 1
    average_time = total_api_time / api_call_count if api_call_count else 0.0

    content_size = len(resp.content) if resp.content else 0

    print(
        # f"[API CALL] {method.upper()} {url} | "
        f"Status: {resp.status_code} | "
        # f"Time: {duration:.4f}s | "
        # f"Size: {content_size} bytes | "
        # f"Cumulative Time: {total_api_time:.4f}s | "
        # f"Calls: {api_call_count} | "
        # f"Avg Time/Call: {average_time:.4f}s"
    )

    # Wait a bit to maintain ~9 requests per second
    # If this sleep doesn't suit your needs, adjust it dynamically or use a token bucket algorithm.
    time.sleep(RATE_LIMIT_INTERVAL)

    return resp