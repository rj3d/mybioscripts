class CachedWriter:
    self.__to_write = []
    self.__cache_cur = 0
    self.__cache_max = 512000000
    self.__cache_filename = ''
    
    def __init__(self, filename):
        self.__cache_filename = filename
        f = open(filename,'w')
        f.close()

    def write(self, string):
        self.__to_write.append(string)
        self.__cache_cur += len(string)
        if self.__cache_cur >= self.__cache_max:
            self.flush()

    def flush(self):
        f = open(self.__cache_filename,'a')
        for line in self.__to_write:
            f.write(line)
        f.close()
        self.__cache_cur = 0
        del self.__to_write[:]
