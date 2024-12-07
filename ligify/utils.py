from collections import defaultdict
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