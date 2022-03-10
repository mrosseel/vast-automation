import sys
from urllib.parse import urlparse
from pathlib import Path
import os


def process_site(url):
    o = urlparse(url)
    p = Path(o.path[1:])
    cwd = Path.cwd()
    varname = p.name
    root = cwd / varname
    postroot = Path(root, p)
    target = postroot.parents[0]
    print(f"o:{o}\np:{p}\np.name:{p.name}\no.netloc{o.netloc}\nroot={root}\nvarname={varname}\npostroot:{postroot}")
    # removing index.html
    print("Removing index.html")
    os.remove(root / "index.html")
    # removing one path level
    # rename /posts/bla/bla to /posts/temp
    os.rename(postroot, postroot.parents[1] / 'temp')
    # remove /posts/bla
    os.rmdir(target)
    # rename posts/temp posts/bla
    os.rename(postroot.parents[1] / 'temp', target)

    # replace ../../.. with ../../
    realindex = Path(target/"index.html")
    txt = realindex.read_text()
    txt = txt.replace(r'../../../', r'../../')
    realindex.write_text(txt)



if __name__ == "__main__":
    process_site(sys.argv[1])
