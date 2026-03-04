#x^2 - dy^2 = 1
#HCC (hyperbolic curve operation)
    #check (check p is on the hcc or not) 
    #findcoord (find the coordinates of giving x on the curve)
    #recordNcoord (print all N coordinates)
    #drawNcoord (draw all N coordinates)
    #add (modadd) 
    #mul (montgomary ladder) 
    #recordncoord (record all n coordinates)
    #drawncoord (draw all n coordinates)
#example
from Crypto.Util.number import isPrime
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator 

class HCC_Discrete_Image:
    def __init__(self, d, p):
        self.zero = (1, 0)
        self.d = d
        self.p = p
        #to check d is perfect square or not if yes curve will degrade
        if pow(d, (p - 1) // 2, p) == 1 or not isPrime(p):
            raise ValueError("The Curve has Defects")

        self.countN, self.containerN = self.recordNcoord()

    def check(self, p):
        if p == self.zero:
            return True

        x, y = p

        if not (0 <= x < self.p and 0 <= y < self.p) or \
        not (x**2 - self.d * y**2 - 1) % self.p == 0:
            raise ValueError("Invalid point")
        return True

    def findcoord(self, x):
        if x < self.p:
            #y^2 = (x^2 - 1) / d
            y = (x**2 - 1) * pow(self.d, -1, self.p) % self.p

            for i in range(0, self.p):
                if i**2 % self.p == y:
                    return((x, i), (x, -i % self.p))
            raise ValueError("Not found")
        else:
            raise ValueError("Out of range")

    def recordNcoord(self):
        #the 1 point means the point at infinity 
        containerN = {self.zero}
        
        #to brute force in order to find the total point 
        for x in range(self.p):
            try:
                p1, p2 = self.findcoord(x)
                #if given x produce 2 same points + 1
                if p1 == p2:
                    containerN.add(p1)
                #if given x produce 2 different points + 2
                else:
                    containerN.add(p1)
                    containerN.add(p2)
            except Exception:
                pass
        return len(containerN), containerN

    def drawNcoord(self):
        x_coord = [p[0] for p in self.containerN if p is not None] 
        y_coord = [p[1] for p in self.containerN if p is not None] 

        #create background
        fig, ax = plt.subplots(figsize=(6, 6))

        #font & title
        title_str = f"$x^2 {-self.d:+}y^2 = 1$ (mod {self.p})"
        ax.set_title(title_str, fontsize=15, pad=15)

        #set the size of background
        margin = self.p * 0.05
        ax.set_xlim(-margin, (self.p - 1) + margin)
        ax.set_ylim(-margin, (self.p - 1) + margin)

        #draw the points
        ax.scatter(x_coord, y_coord, s=60, facecolors='none', edgecolors='#3498db', linewidth=1.5)
        
        #draw the grid
        ax.grid(True, linestyle='-', alpha=0.3)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    #hcc addition
    def add(self, p1, p2):
        self.check(p1)
        self.check(p2)
        #to check p1 / p2 is infinite point or not
        if p1 == self.zero:
            return p2
        elif p2 == self.zero:
            return p1

        x1, y1 = p1
        x2, y2 = p2

        #no need to check inverse it will handle automatically
        x3 = (x1 * x2 + self.d * y1 * y2) % self.p
        y3 = (x1 * y2 + x2 * y1) % self.p

        return (x3, y3)

    #hcc multiplication
    def mul(self, g, n): 
        self.check(g)
        #montgomery ladder
        r0 = self.zero
        r1 = g

        for i in bin(n)[2:]:
            if i == "1":
                r0 = self.add(r0, r1)
                r1 = self.add(r1, r1)
            else:
                r1 = self.add(r1, r0)
                r0 = self.add(r0, r0)
        
        return r0

    def recordncoord(self, g):
        self.check(g)
        containern = []
        current = g

        for i in range(1, self.countN + 1):
            #n will always be a cycle so when n reach infinite point it means n is found
            #to keep letting g adding it self until it reach infinite point
            containern.append(current)
            if current == self.zero:
                return i, containern
            current = self.add(current, g)
        raise ValueError("Invalid order") 

    def drawncoord(self, g):
        countn, containern = self.recordncoord(g)
        x_coord = [p[0] for p in containern if p is not None] 
        y_coord = [p[1] for p in containern if p is not None] 

        #create background
        fig, ax = plt.subplots(figsize=(6, 6))

        #font & title
        title_str = f"Generates Point: {g}"
        ax.set_title(title_str, fontsize=15, pad=15)

        #set the size of background
        margin = self.p * 0.05
        ax.set_xlim(-margin, (self.p - 1) + margin)
        ax.set_ylim(-margin, (self.p - 1) + margin)

        #draw the points
        ax.scatter(x_coord, y_coord, s=60, facecolors='none', edgecolors='#00A877', linewidth=1.5)
        
        #draw the grid
        ax.grid(True, linestyle='-', alpha=0.3)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

#example
if __name__ == "__main__":
    #d = -7
    #p = 19, 97, 131, 419
    d, p = -7, 19
    g = (1, 0)
    hcc = HCC_Discrete_Image(d=d, p=p) 
    print("-------------------------------------------\n") 
    
    countN, containerN = hcc.countN, hcc.containerN
    print(f"HCC: x^2 - {d}y^2 mod({p}) = 1")
    print(f"Counts: {countN}, Container: {containerN}")
    hcc.drawNcoord()
    print("\n-------------------------------------------\n") 

    countn, containern = hcc.recordncoord(g)
    print(f"Generates Point: {g}")
    print(f"Counts: {countn}, Container: {containern}")
    print(f"h: {countN//countn}")
    hcc.drawncoord(g)
    print("\n-------------------------------------------") 
    plt.show()