#!/usr/bin/env python

import sys
import string
import subprocess

#----------------------------------------------------------------------
# Define some handy functions
#----------------------------------------------------------------------
def get_contents(file):
    f = open(file,'r')
    fc = f.readlines()
    f.close()
    return fc

def write_to_file(file,fc):
    f = open(file,'w')
    f.writelines(fc)
    f.close()

#----------------------------------------------------------------------
# The main function
#----------------------------------------------------------------------
def main():
    # Read in identities of dihedrals to change, step size and number of steps per dihedral
    # Note that we will take one extra step per dihedral to capture the initial conformation
    # Read in name of original pdb file, and store the base file name (without the pdb)
    pdb_file = sys.argv[1]
    base = pdb_file.split('/')[-1].split('.')[0]
    gzmat_file = base + ".gzmat"
    diheds = []
    stepsizes = []
    nsteps = []
    for i in range(0,n_dihed):
        diheds.append(sys.argv[3*i+2])
        stepsizes.append(float(sys.argv[3*i+3]))
        nsteps.append(int(sys.argv[3*i+4])+1)

    #----------------------------------------------------------------------
    # Generate gzmat file from pdb file
    #----------------------------------------------------------------------
    process = subprocess.Popen("babel -ipdb {0} -ogzmat {1}".format(pdb_file, gzmat_file).split(), stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
    process.communicate()
    gzmat = get_contents(gzmat_file)
    
    #----------------------------------------------------------------------
    # Find the line numbers and initial values for each dihedral angle in the gzmat file
    #----------------------------------------------------------------------
    lines = []
    values = []
    for i in range(0,n_dihed):
        dihed = diheds[i]
        for j in range(0,len(gzmat)):
            if string.find(gzmat[j],dihed+"=") != -1:
                lines.append(j) 
                value = gzmat[j].split()[1]
                values.append(float(value))
    
    #----------------------------------------------------------------------
    # The main bit of code that actually drives the process 
    # of generating new input files with altered dihedrals
    # Note: the "for x in range(0,y)" is essentially a do loop,
    # using different values of x at each iteration x = 0,1,2,...,y-1
    # Also note that python starts counting at 0 rather than 1
    #----------------------------------------------------------------------
    for i0 in range(0,nsteps[0]):
        # for first dihedral, set new value to initial + step number*step size
        newvalue0 = values[0] + float(i0)*stepsizes[0]
        gzmat[lines[0]] = diheds[0] + "= " + str(newvalue0) + "\n" 
        new_gzmat_file = base + "_" + str(i0).zfill(3) + ".gzmat"
        new_pdb_file = base + "_" + str(i0).zfill(3) + ".pdb"
        if n_dihed == 1:
            # generate file and convert back to pdb format
            write_to_file(new_gzmat_file,gzmat)
            process = subprocess.Popen("babel -igzmat {0} -opdb {1}".format(new_gzmat_file, new_pdb_file).split(), stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.communicate() 
        else:
            for i1 in range(0,nsteps[1]):
                # for second dihedral, set new value to initial + step number*step size
                newvalue1 = values[1] + float(i1)*stepsizes[1]
                gzmat[lines[1]] = diheds[1] + "= " + str(newvalue1) + "\n" 
                new_gzmat_file = base + "_" + str(i0).zfill(3) + "_" + str(i1).zfill(3) + ".gzmat"
                new_pdb_file = base + "_" + str(i0).zfill(3) + "_" + str(i1).zfill(3) + ".pdb"
                if n_dihed == 2:
                    # generate file and convert back to pdb format
                    write_to_file(new_gzmat_file,gzmat)
                    process = subprocess.Popen("babel -igzmat {0} -opdb {1}".format(new_gzmat_file, new_pdb_file).split(), stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
                    process.communicate()
                else:
                    for i2 in range(0,nsteps[2]):
                        # for third dihedral, set new value to initial + step number*step size
                        newvalue2 = values[2] + float(i2)*stepsizes[2]
                        gzmat[lines[2]] = diheds[2] + "= " + str(newvalue2) + "\n" 
                        new_gzmat_file = base + "_" + str(i0).zfill(3) + "_" + str(i1).zfill(3) + "_" + str(i2).zfill(3) + ".gzmat"
                        new_pdb_file = base + "_" + str(i0).zfill(3) + "_" + str(i1).zfill(3) + "_" + str(i2).zfill(3) +  ".pdb"
                        if n_dihed == 3:
                            # generate file and convert back to pdb format
                            write_to_file(new_gzmat_file,gzmat)
                            process = subprocess.Popen("babel -igzmat {0} -opdb {1}".format(new_gzmat_file, new_pdb_file).split(), stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
                            process.communicate()
                        else:
                            print('Error: more than 3 dihedrals, should not be able to get here')
                            sys.exit()
                            
if __name__ == '__main__':
    #----------------------------------------------------------------------
    # Read in arguments from command line
    #----------------------------------------------------------------------
    
    if (len(sys.argv)-2)%3 != 0 or len(sys.argv) == 2:
        print('Usage: scan_dihedral.py <pdb_file_name> <dihedral name 1> <step size 1> <number of steps 1>') 
        print('       ...  ...  ...  ...  <dihedral name N> <step size N> <number of steps N>')
    else:    
        # Calculate the number of dihedrals to change based on the number of arguments supplied
        n_dihed = (len(sys.argv)-2)/3
        if (n_dihed > 3): 
            print('Changing more than 3 dihedrals at once, are you sure?')
            print('If so, you will need to edit the python script to remove the sys.exit() statement')
            print('And write some more do loops in the main part of the code')
            sys.exit()
        main()
