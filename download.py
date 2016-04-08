import urllib.request as urllib
import threading

class Download(threading.Thread):
    def __init__(self, url, filename=None):
        self.url = url
        self.filename = filename
        threading.Thread.__init__(self)
    def get_file(self):
        tries=0
        try:
            urllib.urlretrieve(self.url, self.filename)
            return
        except:
            if tries < 5:
                tries += 1
            else:
                return
    def open_file(self):
        tries = 0
        try:
            x = urllib.urlopen(self.url,timeout=5).read().decode('utf-8')
            return x
        except:
            if tries < 5:
                tries += 1
            else:
                return None