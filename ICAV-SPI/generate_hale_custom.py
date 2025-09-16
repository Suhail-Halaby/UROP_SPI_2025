import h5py as h5
import numpy as np
import os
import sharpy.utils.algebra as algebra

case_name = 'simple_HALE'
route = os.path.dirname(os.path.realpath(__file__)) + '/'



## EXECUTION
flow = ['BeamLoader',
        'AerogridLoader',
        #'NonLinearStatic',
        #'StaticUvlm',
        'StaticTrim',
        #'StaticCoupled',
        'BeamLoads',
        'BeamPlot',
        'AerogridPlot',
        'DynamicCoupled',
        #'Modal',
        #'LinearAssembler',
        #'AsymptoticStability',
        ]


# if free_flight is False, the motion of the centre of the wing is prescribed.
free_flight = True
if not free_flight:
    case_name += '_prescribed'
    amplitude = 0 * np.pi / 180
    period = 3
    case_name += '_amp_' + str(amplitude).replace('.', '') + '_period_' + str(period)




## Data Input



# 1. Flight Conditions  -------------------------------------------------------------------------------------

# Freestream conditions
u_inf = 15
rho = 1.225

alpha = 0.5 * np.pi / 180
beta = 0
roll = 0

gravity = 'on'

# initial control inputs
cs_deflection = -3.8 * np.pi / 180
rudder_static_deflection = 0.0
rudder_step = 0.0 * np.pi / 180
thrust = 0.3 # N, positive is forward

# gust conditions 
gust_intensity = 0.1
gust_length = 1 * u_inf
gust_offset = 0.5 * u_inf 
gust_component = 1  # 'x --> 0 ', 'y --> 1 ', 'z --> 2'

# 2. Numerical Parameters  ----------------------------------------------------------------------------------

n_step = 5
structural_relaxation_factor = 0.15
relaxation_factor = 0.10
tolerance = 1e-4
fsi_tolerance = 1e-3
num_cores = 2

chord_ref = 0.25 # Goland wing reference chord. Used for reduced frequency normalisation
rom_settings = dict()
rom_settings['algorithm'] = 'mimo_rational_arnoldi'  # reduction algorithm
rom_settings['r'] = 6  # Krylov subspace order
frequency_continuous_k = np.array([0.])  # Interpolation point in the complex plane with reduced frequency units
frequency_continuous_w = 2 * u_inf * frequency_continuous_k / chord_ref
rom_settings['frequency'] = frequency_continuous_w

# 3. Geometry  -------------------------------------------------------------------------------------

# Half-spans
span_main = 2.1
span_tail = 0.5
fin_height = 0.4
length_fuselage = 2.05

# Offsets
fuselage_zoffset = 0.01        # with respect to the fuselage centerline
empennage_xoffset = 2.05-0.30  # with respect to the fuselage nose node
wing_xoffset = 0.75            # with respect to the fuselage nose node

# Wing Parameters
chord_main = 0.25
chord_tail = 0.25
chord_fin = 0.30

# Elastic axis offsets
ea_main = 0.228
ea_tail = 0.42
ea_fin = 0.42

# Polyhedral Parameters
subdivision_main = 0.5
lambda_dihedral_outer = 6 * np.pi / 180
lambda_dihedral_inner = 6 * np.pi / 180

span_main1 = (1.0 - subdivision_main) * span_main
span_main2 = subdivision_main * span_main

# motor mounting angle 
epsilon = 1 * np.pi / 180  # positive is upwards


# 4. Beam properties -------------------------------------------------------------------------------------

toggle = False

if toggle == True :
    # Ref. Beam Properties
    ea = 1e7
    ga = 1e5
    gj = 1e4
    eiy = 5e4
    eiz = 5e6

    # Scaling factors
    sigma = 0.1
    sigma_fuselage = 0.1
    sigma_tail = 0.1
else:
    # Wing Beam Properties
    ga = 6.50E+06	
    ea = 9.10E+06	
    gj = 9.86E+04	
    eiy = 7.57E+02	
    eiz = 1.37E+05

    # Fuselage Properties
    ga_f = 5.89E+06	
    ea_f = 8.24E+06	
    gj_f = 6.70E+02	
    eiy_f = 4.69E+02	
    eiz_f = 4.69E+02

    # Empennage Properties
    ga_emp = 1.92E+06	
    ea_emp = 2.69E+06	
    gj_emp = 3.64E+04	
    eiy_emp = 2.30E+02	
    eiz_emp = 5.09E+04


    # Scaling factors
    sigma = 1
    sigma_fuselage = 4
    sigma_tail = 1


# Sectional Mass properties
m_bar_main = 0.716
j_bar_main = 6.94E-03
m_bar_fuselage = 0.318
j_bar_fuselage = 3.62E-05
m_bar_tail = 0.4
j_bar_tail = 1.29E-03


# 5. Lumped Mass Properties -------------------------------------------------------------------------------------

n_lumped_mass = 16

lumped_mass_nodes = np.zeros((n_lumped_mass,), dtype=int)

lumped_mass = np.zeros((n_lumped_mass,))
lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
lumped_mass_position = np.zeros((n_lumped_mass, 3))

lumped_mass[:] = [0.11,0.343,0.051,0.046,0.012,0.15,0.18,0.046,0.436,0.436,0.18,0.045,0.12,0.2,0.1,0.1]
lumped_mass_position[:,0] = [60,36,80,115,225,225,225,335,390,485,445,554,630,790,910,2030]
lumped_mass_position[:,0] = (lumped_mass_position[:,0]/1000)


# 6. Discretisation Parameters  -------------------------------------------------------------------------------------

# chordiwse panels
m = 4

# spanwise elements
n_elem_multiplier = 2
n_elem_main = int(4 * n_elem_multiplier)
n_elem_tail = int(2 * n_elem_multiplier)
n_elem_fin = int(2 * n_elem_multiplier)
n_elem_fuselage = int(2 * n_elem_multiplier)
n_surfaces = 5

# number of nodes per element
n_node_elem = 3

# temporal discretisation
physical_time = 15
tstep_factor = 1.
dt = 1.0 / m / u_inf * tstep_factor
n_tstep = round(physical_time / dt)


##--------------------------------------------------------------------------------------------------------------
## END OF Data Input ##
##--------------------------------------------------------------------------------------------------------------





# 1. Element and Node Processing  -------------------------------------------------------------------------------------

n_elem_main1 = round(n_elem_main * (1 - subdivision_main))
n_elem_main2 = n_elem_main - n_elem_main1

# tallying total number of elements

n_elem = 0
n_elem += 2*n_elem_main1  # inner wings
n_elem += 2*n_elem_main2  # outer wings
n_elem += n_elem_fuselage # fuselage
n_elem += n_elem_fin                # fin
n_elem += 2*n_elem_tail  # tailplane, both halves

# allocating nodes per part

n_node_main1 = n_elem_main1 * (n_node_elem - 1) + 1  # inner wing
n_node_main2 = n_elem_main2 * (n_node_elem - 1) + 1  # outer wing
n_node_main = n_node_main1 + n_node_main2 - 1        # total number of nodes in the main wing

n_node_fuselage = n_elem_fuselage * (n_node_elem - 1) + 1   # fuselage
n_node_fin = n_elem_fin * (n_node_elem - 1) + 1             # fin
n_node_tail = n_elem_tail * (n_node_elem - 1) + 1           # tailplane

# tallying total number of nodes

n_node = 0 
n_node += 2*n_node_main1 - 1        # inner wings
n_node += 2*(n_node_main2 - 1)      # outer wings
n_node += n_node_fuselage - 1       # fuselage
n_node += n_node_fin - 1            # fin
n_node += 2*(n_node_tail - 1)       # tailplane, both halves


# 2. Stiffness and Mass Matrices  -------------------------------------------------------------------------------------

# Stiffnesses


n_stiffness = 3

if toggle == True:
    base_stiffness_main = sigma * np.diag([ea, ga, ga, gj, eiy, eiz])
    base_stiffness_fuselage = base_stiffness_main.copy() * sigma_fuselage
    base_stiffness_tail = base_stiffness_main.copy() * sigma_tail

    base_stiffness_fuselage[4, 4] = base_stiffness_fuselage[5, 5] 
    base_stiffness_tail[4, 4] = base_stiffness_tail[5, 5]
else:
    base_stiffness_main =  sigma*np.diag([ea, ga, ga, gj, eiy, eiz])
    base_stiffness_fuselage = sigma_fuselage*np.diag([ea_f, ga_f, ga_f, gj_f, eiy_f, eiz_f])
    base_stiffness_tail = sigma_tail*np.diag([ea_emp, ga_emp, ga_emp, gj_emp, eiy_emp, eiz_emp])

    base_stiffness_fuselage[4, 4] = base_stiffness_fuselage[5, 5] 
    base_stiffness_tail[4, 4] = base_stiffness_tail[5, 5]

# Masses
n_mass = 3
base_mass_main = np.diag([m_bar_main, m_bar_main, m_bar_main, j_bar_main, 0.5 * j_bar_main, 0.5 * j_bar_main])
base_mass_fuselage = np.diag([m_bar_fuselage,
                              m_bar_fuselage,
                              m_bar_fuselage,
                              j_bar_fuselage,
                              j_bar_fuselage * 0.5,
                              j_bar_fuselage * 0.5])
base_mass_tail = np.diag([m_bar_tail,
                          m_bar_tail,
                          m_bar_tail,
                          j_bar_tail,
                          j_bar_tail * 0.5,
                          j_bar_tail * 0.5])


# 3. Initializing Aero and Structural Properties  -------------------------------------------------------------------------------------

# Beam Geometry
x = np.zeros((n_node,))
y = np.zeros((n_node,))
z = np.zeros((n_node,))

beam_number = np.zeros((n_elem,), dtype=int)
frame_of_reference_delta = np.zeros((n_elem, n_node_elem, 3))
structural_twist = np.zeros((n_elem, 3))
conn = np.zeros((n_elem, n_node_elem), dtype=int)
stiffness = np.zeros((n_stiffness, 6, 6))
elem_stiffness = np.zeros((n_elem,), dtype=int)
mass = np.zeros((n_mass, 6, 6))
elem_mass = np.zeros((n_elem,), dtype=int)
boundary_conditions = np.zeros((n_node,), dtype=int)
app_forces = np.zeros((n_node, 6))

# Aero Geometry
airfoil_distribution = np.zeros((n_elem, n_node_elem), dtype=int)
surface_distribution = np.zeros((n_elem,), dtype=int) - 1
surface_m = np.zeros((n_surfaces,), dtype=int)
m_distribution = 'uniform'
aero_node = np.zeros((n_node,), dtype=bool)
twist = np.zeros((n_elem, n_node_elem))
sweep = np.zeros((n_elem, n_node_elem))
chord = np.zeros((n_elem, n_node_elem,))
elastic_axis = np.zeros((n_elem, n_node_elem,))


#------------------------------------------------------------------------------------------------------------------------
## End of Initialisation, Functions Below ##
#------------------------------------------------------------------------------------------------------------------------


def clean_test_files():
    fem_file_name = route + '/' + case_name + '.fem.h5'
    if os.path.isfile(fem_file_name):
        os.remove(fem_file_name)

    dyn_file_name = route + '/' + case_name + '.dyn.h5'
    if os.path.isfile(dyn_file_name):
        os.remove(dyn_file_name)

    aero_file_name = route + '/' + case_name + '.aero.h5'
    if os.path.isfile(aero_file_name):
        os.remove(aero_file_name)

    solver_file_name = route + '/' + case_name + '.sharpy'
    if os.path.isfile(solver_file_name):
        os.remove(solver_file_name)

    flightcon_file_name = route + '/' + case_name + '.flightcon.txt'
    if os.path.isfile(flightcon_file_name):
        os.remove(flightcon_file_name)

#---------------------------------------------------------------------------------------------------------
## Generate Dyn File ##
#---------------------------------------------------------------------------------------------------------            
            
def generate_dyn_file():
    global dt
    global n_tstep
    global route
    global case_name
    global num_elem
    global num_node_elem
    global num_node
    global amplitude
    global period
    global free_flight

    dynamic_forces_time = None
    with_dynamic_forces = False
    with_forced_vel = False
    if not free_flight:
        with_forced_vel = True

    if with_dynamic_forces:
        f1 = 100
        dynamic_forces = np.zeros((num_node, 6))
        app_node = [int(num_node_main - 1), int(num_node_main)]
        dynamic_forces[app_node, 2] = f1
        force_time = np.zeros((n_tstep,))
        limit = round(0.05 / dt)
        force_time[50:61] = 1

        dynamic_forces_time = np.zeros((n_tstep, num_node, 6))
        for it in range(n_tstep):
            dynamic_forces_time[it, :, :] = force_time[it] * dynamic_forces

    forced_for_vel = None
    if with_forced_vel:
        forced_for_vel = np.zeros((n_tstep, 6))
        forced_for_acc = np.zeros((n_tstep, 6))
        for it in range(n_tstep):
            # if dt*it < period:
            # forced_for_vel[it, 2] = 2*np.pi/period*amplitude*np.sin(2*np.pi*dt*it/period)
            # forced_for_acc[it, 2] = (2*np.pi/period)**2*amplitude*np.cos(2*np.pi*dt*it/period)

            forced_for_vel[it, 3] = 2 * np.pi / period * amplitude * np.sin(2 * np.pi * dt * it / period)
            forced_for_acc[it, 3] = (2 * np.pi / period) ** 2 * amplitude * np.cos(2 * np.pi * dt * it / period)

    if with_dynamic_forces or with_forced_vel:
        with h5.File(route + '/' + case_name + '.dyn.h5', 'a') as h5file:
            if with_dynamic_forces:
                h5file.create_dataset(
                    'dynamic_forces', data=dynamic_forces_time)
            if with_forced_vel:
                h5file.create_dataset(
                    'for_vel', data=forced_for_vel)
                h5file.create_dataset(
                    'for_acc', data=forced_for_acc)
            h5file.create_dataset(
                'num_steps', data=n_tstep)

#---------------------------------------------------------------------------------------------------------
## Generate FEM File ##
#---------------------------------------------------------------------------------------------------------

def generate_fem():
    stiffness[0, ...] = base_stiffness_main
    stiffness[1, ...] = base_stiffness_fuselage
    stiffness[2, ...] = base_stiffness_tail

    mass[0, ...] = base_mass_main
    mass[1, ...] = base_mass_fuselage
    mass[2, ...] = base_mass_tail


    ## Element and Node Arrays
    wn = 0    # written node
    we = 0    # written element
    
    #### NOTE VERIFY BEAM NUMBERING IS CORRECT ####

    ## 1. Fuselage --------------------------------------------------------------------------------------

    # spans nodes 0 to  n_node_fuselage 

    beam_number[we:we + n_elem_fuselage] = 0

    fus_nrecord = wn
    fus_erecord = we

    x[wn:wn + (n_node_fuselage)] = np.linspace(0, length_fuselage, n_node_fuselage)
    z[wn:wn + (n_node_fuselage)] = 0 #fuselage_zoffset

    for ielem in range(n_elem_fuselage):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [0.0, 1.0, 0.0] 
    conn[we, 0] = 0  # fuselage nose node is at the first node 

    elem_stiffness[we:we + n_elem_fuselage] = 1
    elem_mass[we:we + n_elem_fuselage] = 1

    # incrementing counters

    we += n_elem_fuselage           # number of elements written
    wn += n_node_fuselage           # number of nodes written


    # identifying fuselage aft node
    global end_of_fuselage_node
    end_of_fuselage_node = wn - 1

    # identifying wing root node
    fuselage_node_indices = np.arange(fus_nrecord,fus_nrecord+n_node_fuselage)
    fuselage_x = x[fuselage_node_indices]

    global wing_root_node
    wing_root_node = fuselage_node_indices[np.argmin(np.abs(fuselage_x-wing_xoffset))]


    # applied thrust
    # remember this is in B FoR
    app_forces[0] = [-np.cos(epsilon)*thrust, np.sin(epsilon)*thrust, 0, 0, 0, 0]

    
    with open("debug.txt", "w") as f:
        f.write("Fuselage Nodes:    "+str(n_node_fuselage)+"\n")
        f.write("Fuselage Node Positions:\n")
        f.write("Node:"+str(np.arange(fus_nrecord,fus_nrecord+n_node_fuselage))+"\n")
        f.write("x:   "+str(x[0:n_node_fuselage])+ "\n")
        f.write("z:   "+str(z[0:n_node_fuselage])+ "\n")
        f.write("Fuselage Elements: "+str(n_elem_fuselage)+"\n")
        f.write("Connectivities:\n")
        f.write(str(conn[0:n_elem_fuselage, :])+"\n")

        f.write("Wing Root Node: "+str(wing_root_node)+"\n")
        f.write("Wing Root Node Position:\n"+str(x[wing_root_node])+"\n")


    ## 2. Right Wing --------------------------------------------------------------------------------------


    # inner right wing
    # spans n_node_fuselage to n_node_fuselage + n_node_main1


    beam_number[we:we + n_elem_main1] = 1

    rw_nrecord = wn
    rw_erecord = we

    x[wn:wn + n_node_main1-1] = x[wing_root_node]
    y[wn:wn + n_node_main1-1] = np.linspace(0.0, np.cos(lambda_dihedral_inner)*span_main1, n_node_main1)[1:]
    z[wn:wn + n_node_main1-1] = np.linspace(0.0, np.sin(lambda_dihedral_inner)*span_main1, n_node_main1)[1:]
    for ielem in range(n_elem_main1):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]
    conn[we,0] = wing_root_node     
    elem_stiffness[we:we + n_elem_main1] = 0
    elem_mass[we:we + n_elem_main1] = 0
    boundary_conditions[0] = 1

    we += n_elem_main1
    wn += n_node_main1-1

    with open("debug.txt", "a") as f:
        f.write("RW Nodes:"+str(n_node_main1-1)+"\n")
        f.write("RW Node Positions:\n")
        f.write("Node:"+str(np.arange(rw_nrecord,rw_nrecord+n_node_main1-1))+"\n")
        f.write("y:   "+str(y[rw_nrecord:rw_nrecord+n_node_main1-1])+ "\n")
        f.write("z:   "+str(z[rw_nrecord:rw_nrecord+n_node_main1-1])+ "\n")
        f.write("Wing Elements: "+str(n_elem_main1)+"\n")
        f.write("Wing Elements: "+str(n_elem_main)+"\n")
        f.write("Connectivities:\n")
        f.write(str(conn[rw_erecord:rw_erecord+n_elem_main1, :])+"\n")

    rwo_nrecord = wn
    rwo_erecord = we

    # outer right wing
    # spans n_node_fuselage + n_node_main1 to n_node_fuselage + n_node_main1 + n_node_main2 - 1

    beam_number[we:we + n_elem_main1] = 1
    x[wn:wn + n_node_main2 - 1] = x[wing_root_node] 
    y[wn:wn + n_node_main2 - 1] = y[wn - 1] + np.linspace(0.0, np.cos(lambda_dihedral_outer) * span_main2, n_node_main2)[1:]
    z[wn:wn + n_node_main2 - 1] = z[wn - 1] + np.linspace(0.0, np.sin(lambda_dihedral_outer) * span_main2, n_node_main2)[1:]
    for ielem in range(n_elem_main2):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]
    elem_stiffness[we:we + n_elem_main2] = 0
    elem_mass[we:we + n_elem_main2] = 0
    boundary_conditions[wn + n_node_main2 - 2] = -1
    we += n_elem_main2
    wn += n_node_main2 - 1



    with open("debug.txt", "a") as f:
        f.write("RW Nodes:"+str(n_node_main2-1)+"\n")
        f.write("RW Node Positions:\n")
        f.write("Node:"+str(np.arange(rwo_nrecord,rwo_nrecord+n_node_main2-1))+"\n")
        f.write("y:   "+str(y[rwo_nrecord:rwo_nrecord+n_node_main2-1])+ "\n")
        f.write("z:   "+str(z[rwo_nrecord:rwo_nrecord+n_node_main2-1])+ "\n")
        f.write("Wing Elements: "+str(n_elem_main1+n_elem_main2)+"\n")
        f.write("Wing Elements: "+str(n_elem_main)+"\n")
        f.write("Connectivities:\n")
        f.write(str(conn[rwo_erecord:rwo_erecord+n_elem_main2, :])+"\n")



## 3. Left Wing --------------------------------------------------------------------------------------
    
    # inner left wing
    beam_number[we:we + n_elem_main1 - 1] = 2
    
    lw_nrecord = wn
    lw_erecord = we

    x[wn:wn + n_node_main1-1] = x[wing_root_node]
    y[wn:wn + n_node_main1-1] = np.linspace(0.0, -np.cos(lambda_dihedral_inner)*span_main1, n_node_main1)[1:]
    z[wn:wn + n_node_main1-1] = np.linspace(0.0, np.sin(lambda_dihedral_inner)*span_main1, n_node_main1)[1:]
    for ielem in range(n_elem_main1):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [1.0, 0.0, 0.0]
    conn[we, 0] = wing_root_node
    
    elem_stiffness[we:we + n_elem_main1] = 0
    elem_mass[we:we + n_elem_main1] = 0
    we += n_elem_main1
    wn += n_node_main1 - 1


    with open("debug.txt", "a") as f:
        f.write("LW Nodes:"+str(n_node_main1-1)+"\n")
        f.write("LW Node Positions:\n")
        f.write("Node:"+str(np.arange(lw_nrecord,lw_nrecord+n_node_main1-1))+"\n")
        f.write("y:   "+str(y[lw_nrecord:lw_nrecord+n_node_main1-1])+ "\n")
        f.write("z:   "+str(z[lw_nrecord:lw_nrecord+n_node_main1-1])+ "\n")
        f.write("Wing Elements: "+str(n_elem_main1)+"\n")
        f.write("Wing Elements: "+str(n_elem_main)+"\n")
        f.write("Connectivities:\n")
        f.write(str(conn[lw_erecord:lw_erecord+n_elem_main1, :])+"\n")


    # outer left wing
    beam_number[we:we + n_elem_main2] = 2

    lwo_nrecord = wn
    lwo_erecord = we

    x[wn:wn + n_node_main2 - 1] = x[wing_root_node]
    y[wn:wn + n_node_main2 - 1] = y[wn - 1] + np.linspace(0.0, -np.cos(lambda_dihedral_outer) * span_main2, n_node_main2)[1:]
    z[wn:wn + n_node_main2 - 1] = z[wn - 1] + np.linspace(0.0, np.sin(lambda_dihedral_outer) * span_main2, n_node_main2)[1:]
    for ielem in range(n_elem_main2):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [1.0, 0.0, 0.0]
    elem_stiffness[we:we + n_elem_main2] = 0
    elem_mass[we:we + n_elem_main2] = 0
    boundary_conditions[wn + n_node_main2 - 2] = -1
    we += n_elem_main2
    wn += n_node_main2 - 1


    with open("debug.txt", "a") as f:
        f.write("LW Nodes:"+str(n_node_main2-1)+"\n")
        f.write("LW Node Positions:\n")
        f.write("Node:"+str(np.arange(lwo_nrecord,lwo_nrecord+n_node_main2-1))+"\n")
        f.write("y:   "+str(y[lwo_nrecord:lwo_nrecord+n_node_main2-1])+ "\n")
        f.write("z:   "+str(z[lwo_nrecord:lwo_nrecord+n_node_main2-1])+ "\n")
        f.write("Wing Elements: "+str(n_elem_main2)+"\n")
        f.write("Wing Elements: "+str(n_elem_main)+"\n")
        f.write("Connectivities:\n")
        f.write(str(conn[lwo_erecord:lwo_erecord+n_elem_main2, :])+"\n")
        f.write("Connectivity matrix with fuselage and wings")
        f.write(str(conn)+"\n")


# 4. Empennage ----------------------------------------------------------------------------------------------------------

 
   # fin
    beam_number[we:we + n_elem_fin] = 3
    x[wn:wn + n_node_fin - 1] = x[end_of_fuselage_node] #empennage_xoffset
    z[wn:wn + n_node_fin - 1] = z[end_of_fuselage_node] + np.linspace(0.0, fin_height, n_node_fin)[1:]
    for ielem in range(n_elem_fin):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]
    conn[we, 0] = end_of_fuselage_node
    elem_stiffness[we:we + n_elem_fin] = 2
    elem_mass[we:we + n_elem_fin] = 2
    we += n_elem_fin
    wn += n_node_fin - 1
    end_of_fin_node = wn - 1

    with open("debug.txt", "a") as f:
        f.write("Connectivity Matrix With Fin\n")
        f.write(str(conn)+"\n")

    # right tail
    beam_number[we:we + n_elem_tail] = 4
    x[wn:wn + n_node_tail - 1] = x[end_of_fin_node] #empennage_xoffset
    y[wn:wn + n_node_tail - 1] = np.linspace(0.0, span_tail, n_node_tail)[1:]
    z[wn:wn + n_node_tail - 1] = z[end_of_fin_node]
    for ielem in range(n_elem_tail):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]
    conn[we, 0] = end_of_fin_node
    elem_stiffness[we:we + n_elem_tail] = 2
    elem_mass[we:we + n_elem_tail] = 2
    boundary_conditions[wn + n_node_tail - 2] = -1
    we += n_elem_tail
    wn += n_node_tail - 1

    with open("debug.txt", "a") as f:
        f.write("Connectivity Matrix With Right Tailplane\n")
        f.write(str(conn)+"\n")

    # left tail
    beam_number[we:we + n_elem_tail] = 5
    x[wn:wn + n_node_tail - 1] = x[end_of_fin_node] #empennage_xoffset
    y[wn:wn + n_node_tail - 1] = np.linspace(0.0, -span_tail, n_node_tail)[1:]
    z[wn:wn + n_node_tail - 1] = z[end_of_fin_node]
    for ielem in range(n_elem_tail):
        conn[we + ielem, :] = ((np.ones((3,)) * (we + ielem) * (n_node_elem - 1)) +
                               [0, 2, 1])
        for inode in range(n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [1.0, 0.0, 0.0]
    conn[we, 0] = end_of_fin_node
    elem_stiffness[we:we + n_elem_tail] = 2
    elem_mass[we:we + n_elem_tail] = 2
    boundary_conditions[wn + n_node_tail - 3] = -1
    we += n_elem_tail
    wn += n_node_tail - 1

    with open("debug.txt", "a") as f:
        f.write("Connectivity Matrix With Left Tailplane\n")
        f.write(str(conn)+"\n")




# ---------------------------------------------------------------------------------------------------------
## Debugging Scritpt ##
#---------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------

    with h5.File(route + '/' + case_name + '.fem.h5', 'a') as h5file:
        coordinates = h5file.create_dataset('coordinates', data=np.column_stack((x, y, z)))
        conectivities = h5file.create_dataset('connectivities', data=conn)
        num_nodes_elem_handle = h5file.create_dataset(
            'num_node_elem', data=n_node_elem)
        num_nodes_handle = h5file.create_dataset(
            'num_node', data=n_node)
        num_elem_handle = h5file.create_dataset(
            'num_elem', data=n_elem)
        stiffness_db_handle = h5file.create_dataset(
            'stiffness_db', data=stiffness)
        stiffness_handle = h5file.create_dataset(
            'elem_stiffness', data=elem_stiffness)
        mass_db_handle = h5file.create_dataset(
            'mass_db', data=mass)
        mass_handle = h5file.create_dataset(
            'elem_mass', data=elem_mass)
        frame_of_reference_delta_handle = h5file.create_dataset(
            'frame_of_reference_delta', data=frame_of_reference_delta)
        structural_twist_handle = h5file.create_dataset(
            'structural_twist', data=structural_twist)
        bocos_handle = h5file.create_dataset(
            'boundary_conditions', data=boundary_conditions)
        beam_handle = h5file.create_dataset(
            'beam_number', data=beam_number)
        app_forces_handle = h5file.create_dataset(
            'app_forces', data=app_forces)
        lumped_mass_nodes_handle = h5file.create_dataset(
            'lumped_mass_nodes', data=lumped_mass_nodes)
        lumped_mass_handle = h5file.create_dataset(
            'lumped_mass', data=lumped_mass)
        lumped_mass_inertia_handle = h5file.create_dataset(
            'lumped_mass_inertia', data=lumped_mass_inertia)
        lumped_mass_position_handle = h5file.create_dataset(
            'lumped_mass_position', data=lumped_mass_position)
        


#---------------------------------------------------------------------------------------------------------
## Generate Aero File ##
#---------------------------------------------------------------------------------------------------------

def generate_aero_file():

    global x, y, z

# 1. Control Surfaces -------------------------------------------------------------------------------------
    
    # control surfaces
    n_control_surfaces = 2
    control_surface = np.zeros((n_elem, n_node_elem), dtype=int) - 1
    control_surface_type = np.zeros((n_control_surfaces,), dtype=int)
    control_surface_deflection = np.zeros((n_control_surfaces,))
    control_surface_chord = np.zeros((n_control_surfaces,), dtype=int)
    control_surface_hinge_coord = np.zeros((n_control_surfaces,), dtype=float)

    # control surface type 0 = static
    # control surface type 1 = dynamic

    # elevator
    control_surface_type[0] = 0
    control_surface_deflection[0] = cs_deflection
    control_surface_chord[0] = m/4
    control_surface_hinge_coord[0] = 0.25  # nondimensional wrt elastic axis (+ towards the trailing edge)

    # rudder
    control_surface_type[1] = 0
    control_surface_deflection[1] = rudder_static_deflection
    control_surface_chord[1] = 1
    control_surface_hinge_coord[1] = -0.  # nondimensional wrt elastic axis (+ towards the trailing edge)



# 2. Right wing ----------------------------------------------------------------------------------------------
   
   
    we = 0
    wn = 0

    we += n_elem_fuselage
    wn += n_node_fuselage

    #we += n_elem_fuselage           # number of elements written
    #wn += n_node_fuselage           # number of nodes written

    aero_node[wing_root_node] = True  # wing root node is an aero node

    # right wing (surface 0, beam 1)
    i_surf = 0
    airfoil_distribution[we:we + n_elem_main, :] = 0
    surface_distribution[we:we + n_elem_main] = i_surf
    surface_m[i_surf] = m
    aero_node[wn:wn + n_node_main-1] = True
    temp_chord = np.linspace(chord_main, chord_main, n_node_main)
    temp_sweep = np.linspace(0.0, 0 * np.pi / 180, n_node_main)
    node_counter = 0
    for i_elem in range(we, we + n_elem_main):
        for i_local_node in range(n_node_elem):
            if not i_local_node == 0:
                node_counter += 1
            chord[i_elem, i_local_node] = temp_chord[node_counter]
            elastic_axis[i_elem, i_local_node] = ea_main
            sweep[i_elem, i_local_node] = temp_sweep[node_counter]

    with open("debug.txt", "a") as f:
        f.write("\n Aero Surfaces\n")

        f.write("Right Wing Nodes:"+str(n_node_main-1)+"\n")
        f.write("Right Wing Node Indicies:\n")
        f.write("Node:"+str(np.arange(wn, wn+n_node_main-1))+"\n")
        f.write("Right Wing Elements:"+str(n_elem_main)+"\n")
        f.write("Node:"+str(np.arange(we, we+n_elem_main))+"\n")
        f.write("Aero Nodes:\n")
        f.write("Node:"+str(np.arange(wn, wn+n_node_main-1))+"\n")
        f.write(str(aero_node)+"\n")
        f.write("Chord:\n")
        f.write(str(chord)+"\n")
        



    we += n_elem_main
    wn += n_node_main-1




# 3. Left wing ----------------------------------------------------------------------------------------------

    # left wing (surface 1, beam 2)
    i_surf = 1
    airfoil_distribution[we:we + n_elem_main, :] = 0
    # airfoil_distribution[wn:wn + n_node_main - 1] = 0
    surface_distribution[we:we + n_elem_main] = i_surf
    surface_m[i_surf] = m
    aero_node[wn:wn + n_node_main - 1] = True
    # chord[wn:wn + num_node_main - 1] = np.linspace(main_chord, main_tip_chord, num_node_main)[1:]
    # chord[wn:wn + num_node_main - 1] = main_chord
    # elastic_axis[wn:wn + num_node_main - 1] = main_ea
    temp_chord = np.linspace(chord_main, chord_main, n_node_main)
    node_counter = 0
    for i_elem in range(we, we + n_elem_main):
        for i_local_node in range(n_node_elem):
            if not i_local_node == 0:
                node_counter += 1
            chord[i_elem, i_local_node] = temp_chord[node_counter]
            elastic_axis[i_elem, i_local_node] = ea_main
            sweep[i_elem, i_local_node] = -temp_sweep[node_counter]

    with open("debug.txt", "a") as f:

        f.write("Left Wing Nodes:"+str(n_node_main-1)+"\n")
        f.write("Left Wing Node Indicies:\n")
        f.write("Node:"+str(np.arange(wn, wn+n_node_main-1))+"\n")
        f.write("Right Wing Elements:"+str(n_elem_main)+"\n")
        f.write("Node:"+str(np.arange(we, we+n_elem_main))+"\n")
        f.write("Aero Nodes:\n")
        f.write("Node:"+str(np.arange(wn, wn+n_node_main-1))+"\n")
        f.write(str(aero_node)+"\n")
        f.write("Chord:\n")
        f.write(str(chord)+"\n")



    we += n_elem_main
    wn += n_node_main - 1




# 4. Empennage ---------------------------------------------------------------------------------------------

    # # fin (surface 2, beam 3)
    i_surf = 2
    airfoil_distribution[we:we + n_elem_fin, :] = 1
    # airfoil_distribution[wn:wn + n_node_fin] = 0
    surface_distribution[we:we + n_elem_fin] = i_surf
    surface_m[i_surf] = m
    aero_node[end_of_fuselage_node] = True
    aero_node[wn:wn + n_node_fin] = True
    # chord[wn:wn + num_node_fin] = fin_chord
    for i_elem in range(we, we + n_elem_fin):
        for i_local_node in range(n_node_elem):
            chord[i_elem, i_local_node] = chord_fin
            elastic_axis[i_elem, i_local_node] = ea_fin
            control_surface[i_elem, i_local_node] = 1
    # twist[end_of_fuselage_node] = 0
    # twist[wn:] = 0
    # elastic_axis[wn:wn + num_node_main] = fin_ea
    we += n_elem_fin
    wn += n_node_fin - 1
    #
    # # # right tail (surface 3, beam 4)
    i_surf = 3
    airfoil_distribution[we:we + n_elem_tail, :] = 2
    # airfoil_distribution[wn:wn + n_node_tail] = 0
    surface_distribution[we:we + n_elem_tail] = i_surf
    surface_m[i_surf] = m
    # XXX not very elegant
    aero_node[wn:] = True
    # chord[wn:wn + num_node_tail] = tail_chord
    # elastic_axis[wn:wn + num_node_main] = tail_ea
    for i_elem in range(we, we + n_elem_tail):
        for i_local_node in range(n_node_elem):
            twist[i_elem, i_local_node] = -0
    for i_elem in range(we, we + n_elem_tail):
        for i_local_node in range(n_node_elem):
            chord[i_elem, i_local_node] = chord_tail
            elastic_axis[i_elem, i_local_node] = ea_tail
            control_surface[i_elem, i_local_node] = 0

    we += n_elem_tail
    wn += n_node_tail
    #
    # # left tail (surface 4, beam 5)
    i_surf = 4
    airfoil_distribution[we:we + n_elem_tail, :] = 2
    # airfoil_distribution[wn:wn + n_node_tail - 1] = 0
    surface_distribution[we:we + n_elem_tail] = i_surf
    surface_m[i_surf] = m
    aero_node[wn:wn + n_node_tail - 1] = True
    # chord[wn:wn + num_node_tail] = tail_chord
    # elastic_axis[wn:wn + num_node_main] = tail_ea
    # twist[we:we + num_elem_tail] = -tail_twist
    for i_elem in range(we, we + n_elem_tail):
        for i_local_node in range(n_node_elem):
            twist[i_elem, i_local_node] = -0
    for i_elem in range(we, we + n_elem_tail):
        for i_local_node in range(n_node_elem):
            chord[i_elem, i_local_node] = chord_tail
            elastic_axis[i_elem, i_local_node] = ea_tail
            control_surface[i_elem, i_local_node] = 0
    we += n_elem_tail
    wn += n_node_tail

    with h5.File(route + '/' + case_name + '.aero.h5', 'a') as h5file:
        airfoils_group = h5file.create_group('airfoils')
        # add one airfoil
        naca_airfoil_main = airfoils_group.create_dataset('0', data=np.column_stack(
            generate_naca_camber(P=4, M=4)))
        naca_airfoil_tail = airfoils_group.create_dataset('1', data=np.column_stack(
            generate_naca_camber(P=0, M=0)))
        naca_airfoil_fin = airfoils_group.create_dataset('2', data=np.column_stack(
            generate_naca_camber(P=0, M=0)))

        # chord
        chord_input = h5file.create_dataset('chord', data=chord)
        dim_attr = chord_input.attrs['units'] = 'm'

        # twist
        twist_input = h5file.create_dataset('twist', data=twist)
        dim_attr = twist_input.attrs['units'] = 'rad'

        # sweep
        sweep_input = h5file.create_dataset('sweep', data=sweep)
        dim_attr = sweep_input.attrs['units'] = 'rad'

        # airfoil distribution
        airfoil_distribution_input = h5file.create_dataset('airfoil_distribution', data=airfoil_distribution)

        surface_distribution_input = h5file.create_dataset('surface_distribution', data=surface_distribution)
        surface_m_input = h5file.create_dataset('surface_m', data=surface_m)
        m_distribution_input = h5file.create_dataset('m_distribution', data=m_distribution.encode('ascii', 'ignore'))

        aero_node_input = h5file.create_dataset('aero_node', data=aero_node)
        elastic_axis_input = h5file.create_dataset('elastic_axis', data=elastic_axis)

        control_surface_input = h5file.create_dataset('control_surface', data=control_surface)
        control_surface_deflection_input = h5file.create_dataset('control_surface_deflection',
                                                                 data=control_surface_deflection)
        control_surface_chord_input = h5file.create_dataset('control_surface_chord', data=control_surface_chord)
        control_surface_hinge_coord_input = h5file.create_dataset('control_surface_hinge_coord',
                                                                  data=control_surface_hinge_coord)
        control_surface_types_input = h5file.create_dataset('control_surface_type', data=control_surface_type)

#---------------------------------------------------------------------------------------------------------
## Generate NACA Camber Lines ##
#---------------------------------------------------------------------------------------------------------

def generate_naca_camber(M, P):
    mm = M * 1e-2
    p = P * 1e-1

    def naca(x, mm, p):
        if x < 1e-6:
            return 0.0
        elif x < p:
            return mm / (p * p) * (2 * p * x - x * x)
        elif x > p and x < 1 + 1e-6:
            return mm / ((1 - p) * (1 - p)) * (1 - 2 * p + 2 * p * x - x * x)

    x_vec = np.linspace(0, 1, 1000)
    y_vec = np.array([naca(x, mm, p) for x in x_vec])
    return x_vec, y_vec

#---------------------------------------------------------------------------------------------------------
## Generate Solver File ##
#---------------------------------------------------------------------------------------------------------

def generate_solver_file():
    file_name = route + '/' + case_name + '.sharpy'
    settings = dict()
    settings['SHARPy'] = {'case': case_name,
                          'route': route,
                          'flow': flow,
                          'write_screen': 'on',
                          'write_log': 'on',
                          'log_folder': route + '/output/',
                          'log_file': case_name + '.log'}

    settings['BeamLoader'] = {'unsteady': 'on',
                              'orientation': algebra.euler2quat(np.array([roll,
                                                                          alpha,
                                                                          beta]))}
    settings['AerogridLoader'] = {'unsteady': 'on',
                                  'aligned_grid': 'on',
                                  'mstar': int(20 / tstep_factor),
                                  'freestream_dir': [1., 0, 0],
                                  'wake_shape_generator': 'StraightWake',
                                  'wake_shape_generator_input': {'u_inf': u_inf,
                                                                 'u_inf_direction': [1., 0, 0],
                                                                 'dt': dt}}

    settings['NonLinearStatic'] = {'print_info': 'off',
                                   'max_iterations': 150,
                                   'num_load_steps': 1,
                                   'delta_curved': 1e-1,
                                   'min_delta': tolerance,
                                   'gravity_on': gravity,
                                   'gravity': 9.81}

    settings['StaticUvlm'] = {'print_info': 'on',
                              'horseshoe': 'off',
                              'num_cores': num_cores,
                              'n_rollup': 0,
                              'rollup_dt': dt,
                              'rollup_aic_refresh': 1,
                              'rollup_tolerance': 1e-4,
                              'velocity_field_generator': 'SteadyVelocityField',
                              'velocity_field_input': {'u_inf': u_inf,
                                                       'u_inf_direction': [1., 0, 0]},
                              'rho': rho}

    settings['StaticCoupled'] = {'print_info': 'off',
                                 'structural_solver': 'NonLinearStatic',
                                 'structural_solver_settings': settings['NonLinearStatic'],
                                 'aero_solver': 'StaticUvlm',
                                 'aero_solver_settings': settings['StaticUvlm'],
                                 'max_iter': 100,
                                 'n_load_steps': n_step,
                                 'tolerance': fsi_tolerance,
                                 'relaxation_factor': structural_relaxation_factor}

    settings['StaticTrim'] = {'solver': 'StaticCoupled',
                              'solver_settings': settings['StaticCoupled'],
                              'initial_alpha': alpha,
                              'initial_deflection': cs_deflection,
                              'initial_thrust': thrust}

    settings['NonLinearDynamicCoupledStep'] = {'print_info': 'off',
                                               'max_iterations': 950,
                                               'delta_curved': 1e-1,
                                               'min_delta': tolerance,
                                               'newmark_damp': 5e-3,
                                               'gravity_on': gravity,
                                               'gravity': 9.81,
                                               'num_steps': n_tstep,
                                               'dt': dt,
                                               'initial_velocity': u_inf}

    settings['NonLinearDynamicPrescribedStep'] = {'print_info': 'off',
                                                  'max_iterations': 950,
                                                  'delta_curved': 1e-1,
                                                  'min_delta': tolerance,
                                                  'newmark_damp': 5e-3,
                                                  'gravity_on': gravity,
                                                  'gravity': 9.81,
                                                  'num_steps': n_tstep,
                                                  'dt': dt,
                                                  'initial_velocity': u_inf * int(free_flight)}

    relative_motion = 'off'
    if not free_flight:
        relative_motion = 'on'
    settings['StepUvlm'] = {'print_info': 'off',
                            'num_cores': num_cores,
                            'convection_scheme': 2,
                            'gamma_dot_filtering': 6,
                            'velocity_field_generator': 'GustVelocityField',
                            'velocity_field_input': {'u_inf': int(not free_flight) * u_inf,
                                                     'u_inf_direction': [1., 0, 0],
                                                     'gust_shape': '1-cos',
                                                     'gust_parameters': {'gust_length': gust_length,
                                                                         'gust_intensity': gust_intensity * u_inf,
                                                                         'gust_component': gust_component},
                                                     'offset': gust_offset,
                                                     'relative_motion': relative_motion},
                            'rho': rho,
                            'n_time_steps': n_tstep,
                            'dt': dt}

    if free_flight:
        solver = 'NonLinearDynamicCoupledStep'
    else:
        solver = 'NonLinearDynamicPrescribedStep'
    settings['DynamicCoupled'] = {'structural_solver': solver,
                                  'structural_solver_settings': settings[solver],
                                  'aero_solver': 'StepUvlm',
                                  'aero_solver_settings': settings['StepUvlm'],
                                  'fsi_substeps': 200,
                                  'fsi_tolerance': fsi_tolerance,
                                  'relaxation_factor': relaxation_factor,
                                  'minimum_steps': 1,
                                  'relaxation_steps': 150,
                                  'final_relaxation_factor': 0.5,
                                  'n_time_steps': n_tstep,
                                  'dt': dt,
                                  'include_unsteady_force_contribution': 'on',
                                  'postprocessors': ['BeamLoads', 'BeamPlot', 'AerogridPlot'],
                                  'postprocessors_settings': {'BeamLoads': {'csv_output': 'off'},
                                                              'BeamPlot': {'include_rbm': 'on',
                                                                           'include_applied_forces': 'on'},
                                                              'AerogridPlot': {
                                                                  'include_rbm': 'on',
                                                                  'include_applied_forces': 'on',
                                                                  'minus_m_star': 0},
                                                              }}

    settings['BeamLoads'] = {'csv_output': 'off'}

    settings['BeamPlot'] = {'include_rbm': 'on',
                            'include_applied_forces': 'on'}


    settings['AerogridPlot'] = {'include_rbm': 'on',
                                'include_forward_motion': 'off',
                                'include_applied_forces': 'on',
                                'minus_m_star': 0,
                                'u_inf': u_inf,
                                'dt': dt}

    settings['Modal'] = {'print_info': True,
                         'use_undamped_modes': True,
                         'NumLambda': 30,
                         'rigid_body_modes': True,
                         'write_modes_vtk': 'on',
                         'print_matrices': 'on',
                         'continuous_eigenvalues': 'off',
                         'dt': dt,
                         'plot_eigenvalues': False}

    settings['LinearAssembler'] = {'linear_system': 'LinearAeroelastic',
                                   'linear_system_settings': {
                                       'beam_settings': {'modal_projection': False,
                                                         'inout_coords': 'modes',
                                                         'discrete_time': True,
                                                         'newmark_damp': 0.05,
                                                         'discr_method': 'newmark',
                                                         'dt': dt,
                                                         'proj_modes': 'undamped',
                                                         'use_euler': 'on',
                                                         'num_modes': 40,
                                                         'print_info': 'on',
                                                         'gravity': 'on',
                                                         'remove_dofs': []},
                                       'aero_settings': {'dt': dt,
                                                         'integr_order': 2,
                                                         'density': rho,
                                                         'remove_predictor': False,
                                                         'use_sparse': False,
                                                         'remove_inputs': ['u_gust']
                                                         },                                                   
                                        'track_body': 'on'
                                   }}

    settings['AsymptoticStability'] = {'print_info': 'on',
                                       'modes_to_plot': [0,1,2,3,4,5,6],
                                       'display_root_locus': 'on',
                                       'frequency_cutoff': 0,
                                       'export_eigenvalues': 'on',
                                       'num_evals': 100,
                                       #'velocity_analysis': [5,25,20],
                                       }
    




    import configobj
    config = configobj.ConfigObj()
    config.filename = file_name
    for k, v in settings.items():
        config[k] = v
    config.write()


clean_test_files()
generate_fem()
generate_aero_file()
generate_solver_file()
generate_dyn_file()
