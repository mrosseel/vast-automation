from functools import partial
from multiprocessing import Pool

def thing(it, bullshit=False):
    print(f"For {it} we have id {id(it)}")

if __name__ == '__main__':
    func = partial(thing, bullshit=True)
    pool = Pool(4)
    thelist = [x for x in range(1, 20)]
    print([id(x) for x in thelist])
    for entry in pool.imap(func, thelist, chunksize=2):
        pass
