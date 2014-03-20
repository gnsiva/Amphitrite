

class uno():
    # GUI
    def __init__(self,d):
        self.a = 111
        self.dos = d
    def uSetD(self,n):
        self.dos.a = n
    
class dos():
    # Panel
    def __init__(self):
        self.a = 11
        self.tres = tres()
    def check(self):
        print self.a
        
                
class tres():
    # ListCtrl
    def __init__(self):
        self.a = 1
    
    def check(self,ob):
        print ob.a

d = dos()
u = uno(d)
#u.uSetD(5)
d.a = 10

print u.dos.a

