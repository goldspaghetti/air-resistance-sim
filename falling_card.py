from vpython import *
from math import *
#use Lennard-Jones potential and Hashgrids?

dt = 0.01

class Grid():
    def __init__(self, rows=10, columns=10, size=5):
        self.grid = []  # columns, then rows
        self.rows = rows
        self.columns = columns
        self.size = size
        for i in range(columns):
            row = []
            for j in range(rows):
                row.append(ParticleBucket(size=size, pos=vec(i*size, j*size, 0)).populate_bucket())
            self.grid.append(row)
        # self.populate_buckets()
    # def populate_buckets(self):
    #     pass
    def update_buckets(self):   #check all particle positions, first find wrong then change bucket
        to_change = []
        for i in range(self.columns):
            for j in range(self.rows):
                bucket = self.grid[i][j]
                for p in bucket.particles:
                    correct_row = floor(p.s.pos.x/self.size)
                    correct_column = floor(p.s.pos.y/self.size)
                    if correct_row != j or correct_column != i:
                        to_change.append((p, correct_column, correct_row, i, j))
        for info in to_change:
            p = to_change[0]
            self.grid[info[3]][info[4]].pop(p)
            self.grid[info[1]][info[2]].append(p)

class SingleBucket():   # literally a grid with a single bucket (so all particles affect each other)
    def __init__(self, size = 5):
        self.bucket = ParticleBucket(size=size, pos=vec(0, 0, 0))
        # self.bucket.populate_bucket()
    
    def update_buckets(self):
        self.bucket.update_particles()

# class Triangle():
#     def __init__(self):
#         self.bucket = ParticleBucket()

class ParticleBucket():
    def __init__(self, particles=[], size=10, pos=vec(0, 0, 0)):
        self.connections = []   # there should be a better way...
        self.particles = particles
        self.size = size
        self.pos = pos
        # create connections
        self.populate_bucket()
        print(len(self.particles))
        for i in range(len(self.particles)):
            for j in range(len(self.particles)-i-1):
                print(f"{i}: {i}, {len(self.particles)-j-1}")
                self.connections.append(Connection(self.particles[i], self.particles[len(self.particles)-j-1]))   # I hope this works
        print(len(self.connections))

    # def get_particle_connections(self):
    #     for i in range(len(self.particles)):
    #         for j in range(len(self.particles)-i-1):
    #             print(f"{i}: {i}, {len(self.particles)-j-1}")
    #             self.connections.append(Connection(self.particles[i], self.particles[len(self.particles)-j-1]))   # I hope this works
    #     print(len(self.connections))

    def populate_bucket_triangle(self):
        self.particles.append(Particle(pos=vec(0, 0, 0)))
        self.particles.append(Particle(pos=vec(2, 0, 0)))
        self.particles.append(Particle(pos=vec(1, sqrt(3), 0)))

    def populate_bucket(self, density=1):   #particles per 1x1 square (DENSITY NOT IMPLEMENTED YET)
        for i in range(self.size):
            for j in range(self.size):
                self.particles.append(Particle(pos=self.pos+vec(i, j, 0)))

    # updates particle forces, velocity, and position in bucket
    def update_particles(self):
        # for p in self.particles:
        #     p.F = 0
        for c in self.connections:
            c.apply_force()
        for p in self.particles:
            p.update()


class Connection(): # variables shared between two particles
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.dispersion_energy = 1
        self.o = self.p1.radius*2/(2**(1/6))
    
    def apply_force(self):
        r = sqrt((self.p1.s.pos.x-self.p2.s.pos.x)**2 + (self.p1.s.pos.y-self.p2.s.pos.y)**2 + (self.p1.s.pos.z-self.p2.s.pos.z)**2)
        # V = 4*self.dispersion_energy*((self.p1.radius/r)**12-(self.p1.radius/r)**6)
        r_unit = vec(self.p1.s.pos.x-self.p2.s.pos.x, self.p1.s.pos.y-self.p2.s.pos.y, 0).hat
        F = self.dispersion_energy*4*((-12*(self.o/r)**13)+6*((self.o/r)**7))*r_unit
        self.p1.apply_force(-F)
        self.p2.apply_force(F)
        # dV = V-self.past_V
        # d_pos1 = self.p1.s.pos - self.p1.past_pos
        # d_pos2 = self.p2.s.pos - self.p2.past_pos
        # # IGNORING Z AXIS FOR NOW
        # if d_pos1.x == 0:
        #     d_pos1.x = 0.0000001
        # if d_pos1.y == 0:
        #     d_pos1.y = 0.0000001
        # if d_pos2.x == 0:
        #     d_pos2.x = 0.0000001
        # if d_pos2.y == 0:
        #     d_pos2.y = 0.0000001
        # F1 = vec(dV/(d_pos1.x), dV/(d_pos1.y), 0)
        # F2 = vec(dV/d_pos2.x, dV/d_pos2.y, 0)
        # self.p1.apply_force(F1)
        # self.p2.apply_force(F2)
        # self.past_V = V

class Particle():
    def __init__(self, pos=vec(0, 0, 0)):
        self.radius = 0.4
        self.s = sphere(pos=pos, radius=self.radius)
        self.F = vec(0, 0, 0)
        self.v = vec(0.00, 0.000, 0)
        self.m = 1
        # self.past_pos = pos
    # def update_V(self):
    #     self.past_V = 

    def update(self):
        # self.past_pos.x = self.s.pos.x  # to avoid stupid pointer reference
        # self.past_pos.y = self.s.pos.y
        # self.past_pos.z = self.s.pos.z
        a = self.F/self.m
        self.v += dt*a
        self.s.pos += self.v*dt
        self.F = vec(0, 0, 0)
        # self.s.pos += vec(1000, 0, 0)
        # print(f"past: {self.past_pos}, curr: {self.s.pos}")
        #calc new force!

    def apply_force(self, F):
        # print(f"F {F}to {self.s.pos}")
        self.F += F

    # def calc_force(self, p2):
    #     self.past_V = V
    #     r = sqrt((self.s.pos.x-p2.s.pos.x)**2 + (self.s.pos.y-p2.s.pos.y)**2 + (self.s.pos.z-p2.s.pos.z)**2)
    #     V = 4*self.dispersion_energy*((self.radius/r)**12-(self.radius/r)**6)
    #     if self.past_V == None:
    #         self.past_V = V
    #         return 0
    #     dV = V - self.past_V
    #     d_pos = self.pos - 
    #     self.past_V = V
    #     return vec(dV/)

class Card():
    def __init__(self, pos=vec(2.5, 10, 0), l=5, w=1, h=1):
        self.b = box(pos=pos,length=l, width=w, height=h)
        self.v = vec(0, 0, 0)
        self.m=100
        self.F = vec(0, 0, 0)

    def update(self):
        a = self.F/self.m - vec(0, 9.8, 0)
        self.v += a * dt
        self.b.pos += self.v * dt
        self.F = vec(0, 0, 0)

    def process_collision(self, p): # assuming perfectly elastic collision
        # direct calculation of new velocities for both (TOO BORED ::::(((()))))
        # s_dv = 2*p.m*p.v.y/(p.m+self.m)
        # p_dv = (p.m*p.v.y-self.m*self.s_dv)/p.m
        
        # apply a force on both objects (going to assume that the normal force is directly proportional to the displacement (it's probably not but whatever))
        if self.b.pos.x-self.b.length/2 < p.s.pos.x and p.s.pos.x < self.b.pos.x+self.b.length/2 and self.b.height/2 + self.b.pos.y > p.s.pos.y and self.b.pos.y - self.b.height/2 < p.s.pos.y:
            #inside!
            to_right = self.b.length/2 + self.b.pos.x - p.s.pos.x
            to_left = p.s.pos.x + self.b.length/2 - self.b.pos.x
            to_bottom = p.s.pos.y + self.b.height/2 - self.b.pos.y
            to_top = self.b.length/2 + self.b.pos.y - p.s.pos.y
            # print(f"r: {to_right} l: {to_left} b: {to_bottom} t: {to_top}")
            if to_right < to_left and to_right < to_bottom and to_right < to_top:   # I hate this but whatever (takes less time than thinking)
                # force horizontal positive
                F = 10000*to_right
                self.apply_force(vec(F, 0, 0))
                p.apply_force(vec(-F, 0, 0))
            elif to_left < to_right and to_left < to_bottom and to_left < to_top:
                # force horizontal negative
                F = 10000*to_left
                self.apply_force(vec(-F, 0, 0))
                p.apply_force(vec(F, 0, 0))
            elif to_bottom < to_right and to_bottom < to_top and to_bottom < to_left:
                # force vertical negative
                F = 10000*to_bottom
                self.apply_force(vec(0, F, 0))
                p.apply_force(vec(0, -F, 0))
            else:
                # force vertical positive
                F = 10000*to_top
                self.apply_force(vec(0, -F, 0))
                p.apply_force(vec(0, F, 0))
    def apply_force(self, F):
        self.F += F


class Simulation():
    def __init__(self):
        self.t = 0
        self.grid = SingleBucket()
        self.card = Card()


    def run(self):
        while self.t < 50:
            rate(1/dt)
            # print(self.t)
            self.t += dt
            self.grid.update_buckets()
            # for p in self.grid.bucket.particles:
            #     self.card.process_collision(p)
            # self.card.update()

if __name__ == "__main__":
    sim = Simulation()
    sim.run()