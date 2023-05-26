from vpython import *
from math import *
from time import sleep
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
    def __init__(self, length=10, height=20, density=1.5):
        self.bucket = ParticleBucket(length=length, height=height, pos=vec(0, 0, 0), density=density)
        # self.bucket.populate_bucket()
    
    def update_buckets(self):
        self.bucket.update_particles()

# class Triangle():
#     def __init__(self):
#         self.bucket = ParticleBucket()

class ParticleBucket():
    def __init__(self, particles=[], length=20, height=40, density=1.5, pos=vec(0, 0, 0), lennard=False):
        self.connections = []   # there should be a better way...
        self.particles = particles
        # self.size = size
        self.length=length
        self.height=height
        self.pos = pos
        # create connections
        self.populate_bucket(density)
        print(len(self.particles))
        # for i in range(len(self.particles)):
        #     for j in range(len(self.particles)-i-1):
        #         # print(f"{i}: {i}, {len(self.particles)-j-1}")
        #         self.connections.append(Connection(self.particles[i], self.particles[len(self.particles)-j-1]))   # I hope this works
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

    def populate_bucket_test(self):
        self.particles.append(Particle(pos=vec(10, 30, 0), v=vec(0, -50, 0)))
        # self.particles.append(Particle(pos=vec(2, 0, 0), v=vec(-1, 1, 0)))

    def populate_bucket(self, density=1.5):   #particles per 1x1 square (DENSITY NOT IMPLEMENTED YET)
        for i in range(floor(self.length/density)):
            for j in range(floor(self.height/density)):
                self.particles.append(Particle(pos=self.pos+vec(i*density, j*density, 0), density=density))

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
        # if self.past_r == None:
        #     self.past_r = r
        # V = 4*self.dispersion_energy*((self.p1.radius/r)**12-(self.p1.radius/r)**6)
        r_unit = vec(self.p1.s.pos.x-self.p2.s.pos.x, self.p1.s.pos.y-self.p2.s.pos.y, 0).hat
        F = (-1*self.dispersion_energy*4*((-12*((self.o/r)**13))+6*((self.o/r)**7)))*r_unit
        self.p1.apply_force(F)
        self.p2.apply_force(-F)
        # self.past_r = r
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

# assuming 1 atm, giving air around 1.225 kg per m^3, 1.2 gram per liter
# radius is in cm, so radius of 0.5cm is around 0.1 liters, around 0.12g
# volume is in cm^2 * dm (z position is 10 cm now!!!!!), and 0.012 gram per cm^2 * dm
# temp is average kinetic energy, 
class Particle():
    def __init__(self, pos=vec(0, 0, 0), temp=298.15, radius_percent=0.3, density=1.5, v=None):
        self.volume = density   # density is actually volume occupied by each particle...
        self.radius = radius_percent*self.volume/2
        self.s = sphere(pos=pos, radius=self.radius)
        self.F = vec(0, 0, 0)
        self.m=self.volume*0.012
        # self.m = 1.5*0.12
        v_hat = vector.random()
        v_hat.z = 0
        if v_hat.x == 0 and v_hat.y == 0:
            v_hat.x = 1
        self.v = sqrt(temp*1.380/self.m)*v_hat.hat
        if v != None:
            self.v = v
        # print(f"start: {self.v}")
        # self.past_pos = pos
    # def update_V(self):
    #     self.past_V = 

    def update(self):
        # self.past_pos.x = self.s.pos.x  # to avoid stupid pointer reference
        # self.past_pos.y = self.s.pos.y
        # self.past_pos.z = self.s.pos.z
        a = self.F/self.m
        a = vec(0, 0, 0)
        self.v += dt*a
        # print(f"v: {self.v}")
        self.s.pos += self.v*dt
        self.F = vec(0, 0, 0)
        # self.s.pos += vec(1000, 0, 0)
        # print(f"past: {self.past_pos}, curr: {self.s.pos}")
        #calc new force!

    # def calc_partcile_force(self, particles):
        
    #     for p in particles:
    #         r = sqrt((self.s.pos.x-p.s.pos.x)**2 + (self.s.pos.y-p.s.pos.y)**2 + (self.s.pos.z-p.s.pos.z)**2)

    def process_collision(self, particles):
        for p in particles:
            if p != self:
                # distance to other particle
                d = (self.s.pos-p.s.pos).mag
                if d < self.radius*2:
                    # print(f"collision occured between {self.s.pos}, {p.s.pos}")
                    # collided!, assuming mass is same
                    collision_vec = p.s.pos-self.s.pos
                    v1_dir = proj(self.v, collision_vec)
                    v2_dir = proj(p.v, -1*collision_vec)
                    self.v = self.v - v1_dir + v2_dir
                    p.v = p.v - v2_dir + v1_dir
                    p.s.pos = self.s.pos + self.radius*2*(collision_vec.hat)
                    # print(f"c: {self.s.pos}, v: {v2_dir} | {p.s.pos}, v: {v1_dir}")

    # def process_lennard(self, particles):
    #     for p in particles:
    #         r = 

    def apply_force(self, F):
        # print(f"F {F}to {self.s.pos}"
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

class Boundary():   # assumes perfectly elastic collisions with top and bottom, periodic boundary on left and right walls
    def __init__(self, bottom_left=vec(-10, -5, 0), top_right=vec(10, 10, 0)):
        self.bottom_left = bottom_left
        self.top_right = top_right
        self.left_wall = box(pos=vec(bottom_left.x, (bottom_left.y+top_right.y)/2, 0), length=0.5, height=(top_right.y-bottom_left.y))
        self.bottom_wall = box(pos=vec((bottom_left.x+top_right.x)/2, bottom_left.y, 0), length=(top_right.x-bottom_left.x), height=0.5)
        self.top_wall = box(pos=vec((bottom_left.x+top_right.x)/2, top_right.y, 0), length=(top_right.x-bottom_left.x), height=0.5)
        self.right_wall = box(pos=vec(top_right.x, (top_right.y+bottom_left.y)/2, 0), length=0.5, height=(top_right.y-bottom_left.y))

    def process_particle(self, p):
        # perfectly elastic collision with wall
        if p.s.pos.y < self.bottom_left.y:
            if p.v.y < 0:
                p.v.y *= -1
                p.s.pos.y = self.bottom_left.y
        elif p.s.pos.y > self.top_right.y:
            if p.v.y > 0:
                p.v.y *= -1
                p.s.pos.y = self.top_right.y
        if p.s.pos.x > self.top_right.x:
            if p.v.x > 0:
                p.v.x *= -1
                p.s.pos.x = self.top_right.x
        elif p.s.pos.x < self.bottom_left.x:
            if p.v.x < 0:
                p.v.x *= -1
                p.s.pos.x = self.bottom_left.x
        # periodic boundary
        # if p.s.pos.x > self.top_right.x:
        #     p.s.pos.x = self.bottom_left.x+0.01
        # elif p.s.pos.x < self.bottom_left.x:
        #     p.s.pos.x = self.top_right.x-0.01

class Card():
    def __init__(self, pos=vec(2.5, 10, 0), l=5, w=1, h=1, m=8):
        self.b = box(pos=pos,length=l, width=w, height=h)
        self.v = vec(0, 0, 0)
        self.m=m
        self.F = vec(0, 0, 0)

    def update(self):
        a = self.F/self.m - vec(0, 980, 0)
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
            to_top = self.b.height/2 + self.b.pos.y - p.s.pos.y
            # print(f"r: {to_right} l: {to_left} b: {to_bottom} t: {to_top}")
            # if to_top < to_left and to_top < to_bottom and to_top < to_right:
            #     print(f"should collide top")
            if to_right < to_left and to_right < to_bottom and to_right < to_top:   # I hate this but whatever (takes less time than thinking)
                # force horizontal positive
                # print("collide right")
                self.v.x = (2*p.m*p.v.x+self.m*self.v.x-p.m*self.v.x)/(self.m+p.m)
                p.v.x = (2*self.m*self.v.x+p.m*p.v.x-self.m*p.v.x)/(self.m+p.m)
                p.s.pos.x += to_right
                # self.v.x=(2*p.m/(p.m+self.m)*p.v.x-(p.m-self.m)/(p.m+self.m)*self.v.x)
                # p.v.x=(p.m-self.m)/(p.m+self.m)*p.v.x+2*(self.m)/(p.m+self.m)*self.v.x
            elif to_left < to_right and to_left < to_bottom and to_left < to_top:
                # print("collide left")
                # force horizontal negative
                self.v.x = (2*p.m*p.v.x+self.m*self.v.x-p.m*self.v.x)/(self.m+p.m)
                p.v.x = (2*self.m*self.v.x+p.m*p.v.x-self.m*p.v.x)/(self.m+p.m)
                p.s.pos.x -= to_left
                # self.v.x=-1*(2*p.m/(p.m+self.m)*p.v.x-(p.m-self.m)/(p.m+self.m)*self.v.x)
                # p.v.x=-1*((p.m-self.m)/(p.m+self.m)*p.v.x+2*(self.m)/(p.m+self.m)*self.v.x)
            elif to_bottom < to_right and to_bottom < to_top and to_bottom < to_left:
                # force vertical negative
                # print(f"collide bottom")
                self.v.y = (2*p.m*p.v.y+self.m*self.v.y-p.m*self.v.y)/(self.m+p.m)
                p.v.y = (2*self.m*self.v.y+p.m*p.v.y-self.m*p.v.y)/(self.m+p.m)
                p.s.pos.y -= to_bottom
                # self.v.y=-1*(2*p.m/(p.m+self.m)*p.v.y-(p.m-self.m)/(p.m+self.m)*self.v.y)
                # p.v.y=-1*((p.m-self.m)/(p.m+self.m)*p.v.y+2*(self.m)/(p.m+self.m)*self.v.y)
            else:
                # force vertical positive
                # print(f"collide top!")
                # print(f"before: {p.v.y}, {self.v.y}, {p.s.pos.y}")
                self.v.y = (2*p.m*p.v.y+self.m*self.v.y-p.m*self.v.y)/(self.m+p.m)
                p.v.y = (2*self.m*self.v.y+p.m*p.v.y-self.m*p.v.y)/(self.m+p.m)
                p.s.pos.y += to_top
                # print(f"after: {p.v.y}, {self.v.y}, {p.s.pos.y}")
                # self.v.y=2*p.m/(p.m+self.m)*p.v.y-(p.m-self.m)/(p.m+self.m)*self.v.y
                # p.v.y=(p.m-self.m)/(p.m+self.m)*p.v.y+2*(self.m)/(p.m+self.m)*self.v.y
            # if to_right < to_left and to_right < to_bottom and to_right < to_top:   # I hate this but whatever (takes less time than thinking)
            #     # force horizontal positive
            #     F = 10000*to_right
            #     self.apply_force(vec(F, 0, 0))
            #     p.apply_force(vec(-F, 0, 0))
            # elif to_left < to_right and to_left < to_bottom and to_left < to_top:
            #     # force horizontal negative
            #     F = 10000*to_left
            #     self.apply_force(vec(-F, 0, 0))
            #     p.apply_force(vec(F, 0, 0))
            # elif to_bottom < to_right and to_bottom < to_top and to_bottom < to_left:
            #     # force vertical negative
            #     F = 10000*to_bottom
            #     self.apply_force(vec(0, F, 0))
            #     p.apply_force(vec(0, -F, 0))
            # else:
            #     # force vertical positive
            #     F = 10000*to_top
            #     self.apply_force(vec(0, -F, 0))
            #     p.apply_force(vec(0, F, 0))
    def apply_force(self, F):
        self.F += F

    def process_collision_tilted(self, p):
        pass


class Simulation():
    def __init__(self, lennard=False):
        scene.width=1600
        self.t = 0
        self.v_graph = graph(width=400, title="card y velocity vs time", xtitle="time (s)", ytitle="velocity (m/s)", align="left")
        self.g = gcurve(color=color.blue)
        self.g.plot([0, 0])
        self.a_graph = graph(width=400,title="card acceleration vs time", xtitle="time (s)", ytitle="acceleration (m/s^2)", align="left")
        self.a_graph = gcurve(color=color.green)
        self.e_graph = graph(width=400,title="energy of particles vs time", xtitle="time (s)", ytitle="energy (J)", align="left")
        self.e_curve = gcurve(color=color.red)
        self.y_graph = graph(width=400,title="card position vs time", xtitle="time (s)", ytitle="position (m)", align="left")
        self.y_curve = gcurve(color=color.black)
        # self.a = gcurve(color=color.green)
        # scene.width=400
        self.density_test = wtext(text="particle volume")
        self.density_slider = slider(bind=self.set_density, min=1, max=2)
        self.d_val_text = wtext(text="1.5 cm^2*dm")
        scene.append_to_caption('\n')
        self.height_text = wtext(text="card start height(cm)")
        self.height_slider = slider(bind=self.set_height, min=10, max=100)
        self.h_val_text = wtext(text="100 cm")
        scene.append_to_caption('\n')
        self.length_text = wtext(text="card length (cm)")
        self.length_slider = slider(bind=self.set_length, min=1, max=20, step=0.05)
        self.l_val_text = wtext(text="5 cm")
        scene.append_to_caption('\n')
        self.mass_text = wtext(text="card mass (g)")
        self.mass_slider = slider(bind=self.set_mass, min=1, max=50, step=0.1)
        self.m_val_text = wtext(text="3 g")
        scene.append_to_caption('\n')
        self.start_button = button(text="run!", bind=self.run)

        self.length = 5
        self.height = 100
        self.card_mass = 3
        self.ran = False
        self.particle_density = 1.5

        self.completed = False

    def set_mass(self, s):
        self.card_mass = s.value
        self.m_val_text.text = f"{s.value} g"
    def set_length(self, s):
        self.length = s.value
        self.l_val_text.text = f"{s.value} cm"
    def set_height(self, s):
        self.height= s.value
        self.h_val_text.text = f"{s.value} cm"
    def set_density(self, s):
        self.particle_density = s.value
        self.d_val_text.text = f"{s.value} cm^2*dm"

    def run(self):
        self.start_button.disabled=True
        self.grid = SingleBucket(length=40, height=100, density=self.particle_density)
        # self.card = Card(pos=vec(10, self.height, 0), l=self.length, m=self.card_mass)
        self.boundary = Boundary(bottom_left=vec(-1, -1, 0), top_right=vec(40, 100, 0))
        self.card = Card(pos=vec((self.boundary.bottom_left.x+self.boundary.top_right.x)/2, self.height, 0), l=self.length, m=self.card_mass)
        while self.card.b.pos.y >-1:
            rate(1/dt)
            # print(self.t)(
            v_before = self.card.v.mag

            self.t += dt
            for p in self.grid.bucket.particles:
                self.boundary.process_particle(p)
                p.process_collision(self.grid.bucket.particles)
                p.update()
                
            # self.grid.update_buckets()
            for p in self.grid.bucket.particles:
                self.card.process_collision(p)
            self.card.update()
            # adjust card for boundary (shouldn't be here but whatever)
            if self.card.b.pos.x + self.card.b.length/2 > self.boundary.top_right.x and self.card.v.x > 0:
                self.card.v.x *= -1
            if self.card.b.pos.x - self.card.b.length/2 < self.boundary.bottom_left.x and self.card.v.x < 0:
                self.card.v.x *= -1

            E = 0
            for p in self.grid.bucket.particles:
                E += ((p.v.mag)**2)*p.m/2
            # self.g.plot(self.t, E)
            self.g.plot(self.t, self.card.v.y)
            self.y_curve.plot(self.t, self.card.b.pos.y)
            self.a_graph.plot(self.t, (self.card.v.mag-v_before)/dt)
            self.e_curve.plot(self.t, E)
            # self.g.plot(self.t, self.card.b.pos.y)
            # self.a.plot(self.t, self.card.F.mag/self.card.m)
            scene.center = self.card.b.pos
        self.completed = True

    def run_lennard(self):  # not implemented yet, didn't bother to re-implement in the end :(
        self.start_button.disabled=True
        self.grid = SingleBucket(length=40, height=100, density=self.particle_density)
        # self.card = Card(pos=vec(10, self.height, 0), l=self.length, m=self.card_mass)
        self.boundary = Boundary(bottom_left=vec(-1, -1, 0), top_right=vec(40, 100, 0))
        self.card = Card(pos=vec((self.boundary.bottom_left.x+self.boundary.top_right.x)/2, self.height, 0), l=self.length, m=self.card_mass)
        while self.card.b.pos.y >-1:
            rate(1/dt)
            # print(self.t)(
            v_before = self.card.v.mag

            self.t += dt
            for p in self.grid.bucket.particles:
                self.boundary.process_particle(p)
                p.process_collision(self.grid.bucket.particles)
                p.update()
                
            # self.grid.update_buckets()
            for p in self.grid.bucket.particles:
                self.card.process_collision(p)
            self.card.update()
            # adjust card for boundary (shouldn't be here but whatever)
            if self.card.b.pos.x + self.card.b.length/2 > self.boundary.top_right.x and self.card.v.x > 0:
                self.card.v.x *= -1
            if self.card.b.pos.x - self.card.b.length/2 < self.boundary.bottom_left.x and self.card.v.x < 0:
                self.card.v.x *= -1

            E = 0
            for p in self.grid.bucket.particles:
                E += ((p.v.mag)**2)*p.m/2
            # self.g.plot(self.t, E)
            self.g.plot(self.t, self.card.v.y)
            self.y_curve.plot(self.t, self.card.b.pos.y)
            self.a_graph.plot(self.t, (self.card.v.mag-v_before)/dt)
            self.e_curve.plot(self.t, E)
            # self.g.plot(self.t, self.card.b.pos.y)
            # self.a.plot(self.t, self.card.F.mag/self.card.m)
            scene.center = self.card.b.pos
        self.completed = True
            

if __name__ == "__main__":
    sim = Simulation()
    # this sucks, but whatever
    while not sim.completed:
        sleep(0.1)
    # sim.run()