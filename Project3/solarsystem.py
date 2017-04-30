
import numpy as np
get_ipython().magic('matplotlib inline')

import math


class planet:
    
    def __init__(self, mass, position, velocity,name):
        
        self.mass=mass
        self.pos=position
        self.vel=velocity
        self.name=name     
        

    def calc_acc(self, solar_system):
        '''Calculate gravitational acceleration using a=MG/r^2'''
        
        self.acc=np.array([0,0,0])
        G=19.8695*10**(-30)
        
        for plnt in solar_system:#Add up all gravitational acceleration from planets in our solar system
            
            if plnt==self: #Skipping the planet itself
                pass
            
            else:
                
                r= math.sqrt((self.pos[0]-plnt.pos[0])**2+(self.pos[1]-plnt.pos[1])**2)#Find distance between planets
                
                a_mag= G*plnt.mass/(r*r)#Magintude of acceleration
                
                a_dir=(plnt.pos-self.pos)/r #Direction of acceleration
                
                self.acc=self.acc+(a_mag*a_dir)
                
        
class solar_system:
    
    def __init__(self,n_planets=2, max_time=1, step_size=(1/365)):
        
        self.n_planets= n_planets
        self.max_time=max_time
        self.step_size=step_size

        self.system=[]
        
        self.populate()
            
    def populate(self):
        
        '''Values for mass, distance from sun, and tangential velocity taken from https://ssd.jpl.nasa.gov/?planet_phys_par
        and https://nssdc.gsfc.nasa.gov/planetary/factsheet/'''
        
        names=["Sun", "Earth", "Jupiter", "Mercury", "Venus", "Mars", "Saturn", "Uranus", "Neptune", "Pluto"]
        
        masses=[1.989*10**30, 5.97*10**24, 1898*10**24, 0.330*10**24, 4.87*10**24, 0.642*10**24, 568*10**24, 86.8*10**24,102*10**24,0.0146*10**24]
        
        positions=[np.array([0,0,0]), np.array([1,0,0]) ,np.array([5.203,0,0]), np.array([0.39,0,0]), np.array([0.723,0,0]),                    np.array([1.524,0,0]),np.array([9.539,0,0]),np.array([19.18,0,0]),np.array([30.06,0,0]),np.array([39.53,0,0])]
        
        #Proportional to velocity of earth        
        velocities=[np.array([0,0,0]), np.array([0,1,0]), np.array([0,0.434,0]), np.array([0,1.607,0]), np.array([0,1.174,0])                    , np.array([0,0.802,0]), np.array([0,0.323,0]), np.array([0,0.228,0]), np.array([0,0.182,0]), np.array([0,0.159,0])]
        
        
        #Make all our planets. Store them in a list
        for i in range (self.n_planets):
            
            self.system.append(planet(masses[i],positions[i],velocities[i]*6.819943,names[i]))
        
        #Calculate the initial acceleration for each planet based on every other planet in the solar system.
        for plnt in self.system:
            
            plnt.calc_acc(self.system)
            
    def verlet(self):
        
        '''Using velocity verlet and normal verlet equation, find new position, acceleration, and velocity'''
        
        for plnt in self.system:#Update position
            
            plnt.pos=plnt.pos + self.step_size*plnt.vel+(self.step_size**2)/2*plnt.acc
        
        for plnt in self.system:#Update acceleration for new positions
            
            plnt.old_acc=plnt.acc
            plnt.calc_acc(self.system)
        
        for plnt in self.system:#Update velocity for new accelerations
            
            plnt.vel= plnt.vel+self.step_size/2*np.add(plnt.old_acc,plnt.acc)
            
            
    def run(self,p=True ):
        
        
        t=0
        num_runs=0
        
        if p:
            print("Start")

            for plnt in self.system:

                print("{}: Position: {}, Velocity: {}, Acceleration: {}".format(plnt.name, plnt.pos, plnt.vel, plnt.acc))
        
        #While under max_time and haven't exceded number of runs(safeguard)
        while t<self.max_time:
            
            #Increment time by time step
            t+=self.step_size
            num_runs+=1
            
            #Run verlet method
            self.verlet()
            
        if p:

            print("End:")

            for plnt in self.system:

                print("{}: Position: {}, Velocity: {}, Acceleration: {}".format(plnt.name, plnt.pos, plnt.vel, plnt.acc))

            




