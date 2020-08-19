import sys
from urllib.parse import urlparse
from pathlib import Path
import os

def process_site(url):
    o = urlparse(url)
    p = Path(o.path)
    os.rename(o.netloc, p.name)
    os.remove(Path(p.name) / 'posts' / 'index.html')


if __name__ == "__main__":
    process_site(sys.argv[1])

