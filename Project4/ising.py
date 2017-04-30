import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
import random

class ising_model:
    '''A class to run the ising model. Takes in n for nxn square matrix, max number of Monte Carlo cycles per
    temperature, initial and final temperature of the run, and the temperature steps inbetween them '''
    def __init__(self, size=2, max_cycles=20, init_temp=1, final_temp=4, temp_step=0.5 ):
        
        self.size=size
        self.matrix=np.ones((size,size))
        self.mcs=max_cycles
        self.init_temp=init_temp
        self.final_temp=final_temp
        self.temp_step=temp_step
        
        self.E=0
        self.M=0
        
        self.results=[]
        
    def main(self):
        self.temp=self.init_temp
        
        while self.temp<=self.final_temp :
            
            self.E=0
            self.M=0
            
            self.initialize()
            
            average=np.zeros(5)
            
            for cycles in range(self.mcs):
                
                self.metropolis()
                
                average[0]+=self.E; average[1]+=self.E**2
                average[2]+=self.M; average[3]+=self.M**2; average[4]+= np.fabs(self.M)
            
            self.output(average)
            
            self.temp+=self.temp_step
            
            
        self.plot()
        
    def initialize(self):
        '''Sets all matrix elements to 1 below crit temp. Sums up matrix elements to find M. Uses deltaE to find E'''
        for i in range(len(self.matrix)):
            for k in range(len(self.matrix[i])):
                if self.temp < 1.5:
                    self.matrix[i,k]=1
                self.M+=self.matrix[i,k]
                
        for x in range(len(self.matrix)):
            for y in range(len(self.matrix[x])):
                
                self.E-=self.initE(x,y)
                

                
    def initE(self,x,y):                
                dE=0
                
                if x==0:
                    dE+=self.matrix[-1,y]
                else:
                    dE+=self.matrix[x-1,y]
                    
                if y==0:
                    dE+=self.matrix[x,-1]
                else:
                    dE+=self.matrix[x,y-1]
                return self.matrix[x,y]*dE

    def deltaE(self,x,y):
        '''Sums up the spins around given point. Looks at the other end of the board for points on the boarder'''
        dE=0       
        
        if x==0:
            dE+=self.matrix[-1,y]
        else:
            dE+=self.matrix[x-1,y]
            
        if x==self.size-1:
            dE+=self.matrix[0,y]
        else:
            dE+=self.matrix[x+1,y]
                
        if y==0:
            dE+=self.matrix[x,-1]
        else:
            dE+=self.matrix[x,y-1]
            
        if y==self.size-1:
            dE+=self.matrix[x,0]
        else:
            dE+=self.matrix[x,y+1]
            
        return 2*self.matrix[x,y]*dE
            
    def metropolis(self):
        
        for i in range(len(self.matrix)):
            for k in range(len(self.matrix[0])):
                
                x=random.randrange(self.size)
                y=random.randrange(self.size)
                
                dE=self.deltaE(x,y)
                
                r=random.random()
                w=np.exp(-dE/self.temp)
                
                if r<=w:
                    
                    self.matrix[x,y] *= -1
                    
                    self.M+=2*self.matrix[x,y]
                    self.E+=dE
                
    def output(self,average,mtype="temp"):
        
        norm=average/self.mcs
        
        Evariance=(norm[1]-norm[0]**2)/self.size**2
        Mvariance=(norm[3]-norm[2]**2)/self.size**2
        M2variance=(norm[3]-norm[4]**2)/self.size**2
        
        placeholder=[0,1,2,3,4]
        
        if mtype=="temp":
            placeholder[0]=self.temp
        else:
            placeholder[0]=mtype
        placeholder[1]=norm[0]/self.size**2
        placeholder[2]=Evariance/self.temp**2
        placeholder[3]=M2variance/self.temp
        placeholder[4]=norm[4]/self.size**2

        
        self.results.append(placeholder)
        
    def plot(self):
        
        plot_matrix=np.asarray(self.results)
        
        plot_matrix=plot_matrix.transpose()
        
        plot_matrix=np.around(plot_matrix,8)
                
        plt.plot(plot_matrix[0],plot_matrix[1])
        plt.title("Temperature vs Energy")
        plt.xlabel("Temperature(K)")
        plt.ylabel("Energy(J)")
        plt.show()
        
        plt.plot(plot_matrix[0],plot_matrix[2])
        plt.title("Temperature vs Specific Heat")
        plt.xlabel("Temperature(K)")
        plt.ylabel("Specific Heat")
        plt.show()
        
        plt.plot(plot_matrix[0],plot_matrix[3])
        plt.title("Temperature vs Magnetic Susceptibility")
        plt.xlabel("Temperature(K)")
        plt.ylabel("Susceptibility")
        plt.show()
        
        plt.plot(plot_matrix[0],plot_matrix[4])
        plt.title("Temperature vs Magnetic Moment")
        plt.xlabel("Temperature(K)")
        plt.ylabel("Magnetic Moment")
        plt.show()
        
