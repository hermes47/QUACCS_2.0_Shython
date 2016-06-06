""" This file contains a number of functions that may be useful for completing
possible python extension exercises. Some functions are complete, others will
require completion before they will work. To use, copy and paste the entire
function into your script file.
"""
import numpy as np
#----------------------------------------------------------------------
# Given a number of lists, this function will return all combinations
# of those lists. eg: a = [1,2]; b = [3,4]; c = [5,6]
# returns: [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,3,6),(2,4,5),(2,4,6)]
#----------------------------------------------------------------------
def all_combinations(*args):
    import itertools
    return list(itertools.product(*args))

#----------------------------------------------------------------------
# An example of how to obtain output from a command line application
# within the python environment. Prints the length of each line, followed
# by the line itself.
#----------------------------------------------------------------------
def extract_results():
    import subprocess
    process = subprocess.Popen('ls'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    print "Output:"
    for line in output.split('\n'):
        print len(line), line
    print '\nError:'
    for line in error:
        print len(line), line


#----------------------------------------------------------------------
# The following functions provide examples on generating plots using python
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Plot a 1 dimensional graph, loading data from example data
#----------------------------------------------------------------------
def plot_single_dihedral():
    import matplotlib.pyplot as plt
    # Load the data file
    angles, energies = [], []
    with open('../data/1D_energies.txt','r') as f:
        for line in f:
            data = line.split()
            angles.append(float(data[0]))
            energies.append(float(data[1]))
            
    # Generate the plot
    plt.plot(angles, energies, 'rx')
    plt.xlabel("Angle (degrees)")
    plt.ylabel(r"Energy (kJ mol$^{-1}$)")
    plt.title("Dihedral angle rotation energy profile")
    plt.show()
    plt.close()
    
#----------------------------------------------------------------------
# Plot 2D data as a colour contour plot
#----------------------------------------------------------------------
def plot_double_dihedral_flat():
    import matplotlib.pyplot as plt
    import scipy.interpolate
    import numpy as np
    # Load the data file
    phi, psi, energies = [], [], []
    with open('../data/2D_energies.txt','r') as f:
        for line in f:
            data = line.split()
            phi.append(float(data[0]))
            psi.append(float(data[1]))
            energies.append(float(data[2]))

    # Set up a regular grid of interpolation points
    x = np.asarray(phi)
    y = np.asarray(psi)
    z = np.asarray(energies)
    z_min = z.min()
    z -= z_min
    xi, yi = np.linspace(-180, 180, 100), np.linspace(-180, 180, 100)
    xi, yi = np.meshgrid(xi,yi)
    
    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
    zi = rbf(xi,yi)
    
    # Generate the plot
    # comment out the line below and use the line underneath it to show why interpolation is used
    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[-180,180,-180,180])
    #plt.imshow([x, y], vmin=z.min(), vmax=z.max(), origin='lower', extent=[-180,180,-180,180])
    plt.colorbar()
    plt.xlim((-180,180))
    plt.ylim((-180,180))
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\psi$')
    plt.title(r'$\phi$ - $\psi$ plot for dialanine peptide')
    plt.show()
    plt.close()
  
#----------------------------------------------------------------------
# Plot 2D data as a wireframe surface plot
#----------------------------------------------------------------------  
def plot_double_dihedral_surface():
    from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
    import matplotlib.pyplot as plt
    import scipy.interpolate
    import numpy as np
    # Load the data file
    phi, psi, energies = [], [], []
    with open('../data/2D_energies.txt','r') as f:
        for line in f:
            data = line.split()
            phi.append(float(data[0]))
            psi.append(float(data[1]))
            energies.append(float(data[2]))

    # Set up a regular grid of interpolation points
    x = np.asarray(phi)
    y = np.asarray(psi)
    z = np.asarray(energies)
    z_min = z.min()
    z -= z_min
    xi, yi = np.linspace(-180, 180, 100), np.linspace(-180, 180, 100)
    xi, yi = np.meshgrid(xi,yi)
    
    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='cubic')
    zi = rbf(xi,yi)
    # Generate the plot
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot_wireframe(xi, yi, zi, rstride=5, cstride=5)
    plt.xlim((-180,180))
    plt.ylim((-180,180))
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\psi$')
    ax.set_zlabel(r"Energy (kJ mol$^{-1}$)")
    plt.title(r'$\phi$ - $\psi$ surface plot for dialanine peptide')
    plt.show()
    plt.close()


#----------------------------------------------------------------------
# Perform a linear least squares fit to some data
#---------------------------------------------------------------------- 
def fit_torsion_terms_lls(energies, angles):
    #---------------------------------------------------------------------- 
    # To perform linear least squares fitting, one must solve the matrix
    # equation Ax = b for x. This is basically a series of linear equations
    # which can be solved simultaneously. A contains the coefficients, and
    # b the vector of sums. As performing a least squares, the system can be
    # either under, over or exactly defined.
    #---------------------------------------------------------------------- 
    
    assert len(energies) == len(angles)
    m = 6; t = 1
    coefficient_matrix = np.zeros((1,m*t))
    target_vector = np.array((0))
    num_of_coefficients = 6
    # for each data point
    for i in range(len(energies)):
        # build each line of the coefficient matrix
        energy_vector = np.zeros((1,num_of_coefficients))
        
        # add each coefficient to the matrix. Determine what the coefficient
        # is and then place it at the right point in the vector
        for j in range(num_of_coefficients):
            # Calculate the value of the coefficient by some means
            #---------------------------------------------------------------------- 
            # Add your code here
            #----------------------------------------------------------------------
            
            # add to the correct position in the vector 
            energy_vector[0, j] = coefficient_value
        
        # Add the vector of coefficients and the value to the respective array
        coefficient_matrix = np.vstack((coefficient_matrix, energy_vector))
        target_vector = np.hstack((target_vector,np.array((energies[i]))))
        
    # Get rid of the first list of zeros from each, beacuse we don't need them
    coefficient_matrix = coefficient_matrix[1:]
    target_vector = target_vector[1:]
    
    # Fit using QR decomposition
    q,r = np.linalg.qr(coefficient_matrix)
    d = np.dot(np.transpose(q),target_vector)
    x = np.linalg.solve(r, d)
    # Or fit using SVD
    #x = np.linalg.lstsq(matrix, target)[0]
   
   
    # x gives the least squares fit
    return x


def main():
    plot_double_dihedral_surface()
    
if __name__ == "__main__":
    main()
