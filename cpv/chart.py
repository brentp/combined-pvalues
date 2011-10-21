GCE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-.'

def encode(numbers, GCE=GCE, url=True):
    """
    do extended encoding on a list of numbers for the google chart api
    >>> encode([1690, 90,1000])
    'chd=e:aaBaPo'
    """
    encoded = []
    if url:
        # we know they came in scaled between 0 and 1.
        numbers = [n * 100 for n in numbers]

    numbers = [int(n) for n in numbers]

    for number in numbers:
        if number > 4095: raise ValueError('too large')
        first, second = divmod(number, len(GCE))
        encoded.append("%s%s" % (GCE[first], GCE[second]))
    if url:
        return ("chd=e:%s" % ''.join(encoded)), numbers
    else:
        return encoded


def chart(acfs, base="https://chart.googleapis.com/chart?cht=bvs&chs=550x250",
    color="224499"):
    values = [float(acf) for kbin, acf in acfs]
    url, numbers = encode(values)
    url += "&chco=" + color
    xlabels = "|".join("%s-%s" % k for k, v in acfs)
    # draw the x-axis ... and the bin-names
    enc = encode([0, max(numbers)], url=False)
    url += "&chxt=x,y&chxl=0:|" + xlabels

    url += "&chxr=1,0,1|"
    url += "&chm=N,000000,0,-1,12"
    url += "&chbh=40,17,17"

    return base + "&" + url
