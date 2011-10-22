
def encode(numbers):
    """
    do extended encoding on a list of numbers for the google chart api
    >>> encode([1690, 90,1000])
    'chd=e:aaBaPo'
    """
    return "&chd=t:%s" % ','.join(["%.3f" % n for n in numbers])


def chart(acfs, base="https://chart.googleapis.com/chart?cht=bvs&chs=550x250",
    color="224499"):
    values = [float(acf) for kbin, acf in acfs]
    url = encode(values)
    url += "&chco=" + color
    xlabels = "|".join("%s-%s" % k for k, v in acfs)
    # draw the x-axis ... and the bin-names
    url += "&chxt=x,y&chxl=0:|" + xlabels

    url += "&chxr=1,0,%.2f" % (max(values) + 0.02)
    url += "&chm=N,000000,0,-1,12"
    url += "&chbh=40,17,17&chds=a"
    return base + url
