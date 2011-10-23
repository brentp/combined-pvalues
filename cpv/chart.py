
def encode(numbers):
    """
    do extended encoding on a list of numbers for the google chart api
    >>> encode([1690, 90,1000])
    'chd=e:aaBaPo'
    """
    return "&chd=t:%s" % ','.join(["%.3f" % n for n in numbers])


def chart(values, xlabels, base="https://chart.googleapis.com/chart?cht=bvs&",
    color="224499"):
    url = encode(values)
    chart_width = min(1000, 24 + 50 * len(values))
    url += "&chs=%ix250" % chart_width
    url += "&chco=" + color
    # draw the x-axis ... and the bin-names
    url += "&chxt=x,y&chxl=0:|" + xlabels

    url += "&chxr=1,0,%.2f" % (max(values) + 0.02)
    url += "&chm=N,000000,0,-1,12"
    url += "&chbh=42,6,12&chds=a"
    return base + url
