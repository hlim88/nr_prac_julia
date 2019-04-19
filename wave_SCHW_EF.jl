# Solving wave equation on Schwarzschil background in
# the Eddington-Finkelstein coordinates

using PyCall

np = pyimport("numpy")

# Parameters
N = 5
dx = 1/N
dt = 0.005
#cfl = dt/dx
cfl = 1
t_steps = 250

# Define grid
x_grid = np.linspace(0,1,N)

#initialize arrays
phi = np.zeros((t_steps, N))
r = np.zeros((t_steps, N))
s = np.zeros((t_steps, N))

#define initial data
x_c = 0.5
width = 0.05
amp = 0.5
idsignum = 1. # -1=ingoing; 0=t-symmetric; +1=outgoing

r[1, :] = amp * np.exp(-(x_grid.-x_c).^2/width^2)
s[1, :] = idsignum * amp * np.exp(-(x_grid.-x_c).^2/width^2)

# Define matrix A for our discretization scheme

A = np.zeros((2*N, 2*N))
B = np.zeros((2*N, 2*N))

for i = 1:N
    if i == 1
        for j = 1:2*N
            if j==1
                A[i,j] = 1
            end
            if j==N+1
                A[i,j] = cfl/4
            end
            if j==N+2
                A[i,j] = -cfl/4
            end
        end
    elseif i == N
        for j = 1:2*N
            if j==N
                A[i,j] = 1
            end
            if j==2*N-1
                A[i,j] = cfl/4
            end
            if j==2*N
                A[i,j] = -cfl/4
            end
        end
    else 
        for j = 1:2*N
            if j==i
                A[i,j] = 1
            end
            if j==N+i
                A[i,j] = cfl/4
            end
            if j==N+i+2
                A[i,j] = -cfl/4
            end
        end
    end
end
    
for i = N+1:2*N
    if i == N+1
        for j = 1:2*N
            if j==N+1
                A[i,j] = 1
            end
            if j==1
                A[i,j] = cfl/4
            end
            if j==2
                A[i,j] = -cfl/4
            end
        end
    elseif i == 2*N
        for j = 1:2*N
            if j==2*N
                A[i,j] = 1
            end
            if j==N-1
                A[i,j] = cfl/4
            end
            if j==N
                A[i,j] = -cfl/4
            end
        end
    else 
        for j = 1:2*N
            if j==i
                A[i,j] = 1
            end
            if j==i-N
                A[i,j] = cfl/4
            end
            if j==i+2-N
                A[i,j] = -cfl/4
            end
        end
    end
end
    show(A)
    show(A[1,1])

