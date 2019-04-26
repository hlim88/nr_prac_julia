# Solving wave equation on Schwarzschild background in
# the Eddington-Finkelstein coordinates

using PyCall

np = pyimport("numpy")

# Parameters
N = 10
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

# Define matrix B
B = np.zeros((2*N, 2*N))

for i = 1:N
    if i == 1
        for j = 1:2*N
            if j==1
                B[i,j] = 1
            end
            if j==N+1
                B[i,j] = -cfl/4
            end
            if j==N+2
                B[i,j] = cfl/4
            end
        end
    elseif i == N
        for j = 1:2*N
            if j==N
                B[i,j] = 1
            end
            if j==2*N-1
                B[i,j] = -cfl/4
            end
            if j==2*N
                B[i,j] = cfl/4
            end
        end
    else 
        for j = 1:2*N
            if j==i
                B[i,j] = 1
            end
            if j==N+i
                B[i,j] = -cfl/4
            end
            if j==N+i+2
                B[i,j] = cfl/4
            end
        end
    end
end
    
for i = N+1:2*N
    if i == N+1
        for j = 1:2*N
            if j==N+1
                B[i,j] = 1
            end
            if j==1
                B[i,j] = -cfl/4
            end
            if j==2
                B[i,j] = cfl/4
            end
        end
    elseif i == 2*N
        for j = 1:2*N
            if j==2*N
                B[i,j] = 1
            end
            if j==N-1
                B[i,j] = -cfl/4
            end
            if j==N
                B[i,j] = cfl/4
            end
        end
    else 
        for j = 1:2*N
            if j==i
                B[i,j] = 1
            end
            if j==i-N
                B[i,j] = -cfl/4
            end
            if j==i+2-N
                B[i,j] = cfl/4
            end
        end
    end
end


# Update routine
function update(t_steps)
    u = np.zeros(2*N)
      for i = 1:2*N
          if i<N+1
              u[i] = r[t_steps,i]
          else
              u[i] = s[t_steps,i-N]
          end
      end
    return u
end

# Update r and s using the solved-for vector u at t_step
function update_rs(u_sol, t_steps)
    for i = 1:2*N
        if i<N+1
            r[t_steps, i] = u_sol[i]
        else
            s[t_steps, i-N+1] = u_sol[i+1]
            # Applying Sommerfeld condition for boundray
            s[t_steps, 1] = r[t_steps,1]
            s[t_steps, N] = -r[t_steps,N]
        end
    end
end

# Solve the equations
for i = 1:t_steps
    u = update(i)
    rhs = B.*u
    u_sol = np.linalg.solve(A,rhs)
    update_rs(u_sol,i)
end

# Save data
#for i=1:t_steps
#    np.savetxt('result.dat', [x_grid, r[i,:]])
#end
