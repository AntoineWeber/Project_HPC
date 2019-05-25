

#include "kernels.cuh"


__global__ void initialize_particles_uni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICULES)
    {
        mass[i + offset] = (double)MASS;

        x_pos[i + offset] = curand_uniform(&state)*BOUNDS*2 - BOUNDS;
        y_pos[i + offset] = curand_uniform(&state)*BOUNDS*2 - BOUNDS;

        // set velocity to 0
        x_vel[i + offset] = 0.0;
        y_vel[i + offset] = 0.0;

        // set acceleration to 0
        x_acc[i + offset] = 0.0;
        y_acc[i + offset] = 0.0;

        offset += stride;
    }
}

__global__ void initialize_particles_circle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICULES)
    {
        mass[i + offset] = (double)MASS;

        double r = CIRCLE_OFFSET + curand_uniform(&state)*2*BOUNDS;
        double alpha = curand_uniform(&state)*2*M_PI;

        x_pos[i + offset] = r*cos(alpha);
        y_pos[i + offset] = r*sin(alpha);

        // set velocity to 0
        x_vel[i + offset] = 0.0;
        y_vel[i + offset] = 0.0;

        // set acceleration to 0
        x_acc[i + offset] = 0.0;
        y_acc[i + offset] = 0.0;

        offset += stride;
    }

}


__global__ void compute_displacements(Node *d_tree, double* d_x, double* d_y, double* d_vx, double* d_vy, double* d_ax, double* d_ay, double* d_mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    double fx;
    double fy;

    int my_index;

    while (i + offset < N_PARTICULES)
    {
        my_index = i+offset;
        fx = 0.0;
        fy = 0.0;

        // if particle in the far space
        if (d_x[my_index] > FAR_SPACE || d_x[my_index] < -FAR_SPACE || d_y[my_index] > FAR_SPACE || d_y[my_index] < -FAR_SPACE)
        {
            fx = 0.0;
            fy = 0.0;
        }

        else
        {
            /*
            double r;
            for (int j=21; j<85; j++)
            {
                r = sqrt((d_tree[j].m_x_center-d_x[my_index])*(d_tree[j].m_x_center-d_x[my_index]) + (d_tree[j].m_y_center-d_y[my_index])*(d_tree[j].m_y_center-d_y[my_index]));
                r = r>EPSILON ? r : EPSILON;
                fx += (G*d_mass[my_index]*d_tree[j].m_av_mass*(d_tree[j].m_x_center-d_x[my_index]))/r;
                fy += (G*d_mass[my_index]*d_tree[j].m_av_mass*(d_tree[j].m_y_center-d_y[my_index]))/r;
            }
            */
            compute_branch(1, 0, 4, d_tree, d_x[my_index], d_y[my_index], d_mass[my_index], fx, fy);
        }
        
        offset += stride;
    }
    cudaDeviceSynchronize();
    
    offset = 0;
    while (i + offset < N_PARTICULES)
    {
        my_index = i+offset;
        d_ax[my_index] = fx / d_mass[my_index];
        d_ay[my_index] = fy / d_mass[my_index];

        d_vx[my_index] += d_ax[my_index]*TIMESTEP;
        d_vy[my_index] += d_ay[my_index]*TIMESTEP;

        d_x[my_index] += 0.5*d_ax[my_index]*TIMESTEP*TIMESTEP + d_vx[my_index]*TIMESTEP;
        d_y[my_index] += 0.5*d_ay[my_index]*TIMESTEP*TIMESTEP + d_vy[my_index]*TIMESTEP;

        offset += stride;
    }
}

__device__ void compute_branch(int absOff, int depthOff, int nnode, Node *tree, double x, double y, double m, double &fx, double &fy)
{
    double r, s;
    for (unsigned int i=0; i<4; i++)
    {
        // if child exists
        if (tree[absOff + depthOff + i].m_av_mass != 0.0)
        {
            r = sqrt((tree[absOff + depthOff + i].m_x_center-x)*(tree[absOff + depthOff + i].m_x_center-x) \
                                 + (tree[absOff + depthOff + i].m_y_center-y)*(tree[absOff + depthOff + i].m_y_center-y));

            // nan check
            if (!(r==r))
            {
                continue;
            }
            // if external node, compute force with it and add it to force component
            if (!tree[absOff + depthOff + i].hasChildren)
            {
                if (r > EPSILON)
                {
                    fx += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_x_center-x))/r;
                    fy += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_y_center-y))/r;
                }
                else
                {
                    fx += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_x_center-x))/EPSILON;
                    fy += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_y_center-y))/EPSILON;
                }
                
            }
            // if internal node
            else
            {
                // compute the quotient s/r
                s = tree[absOff + depthOff + i].m_s;

                // if too far, consider as a single body
                //std::cout << s/r << " " << std::endl;
                if (s / r < THETA)
                {
                    //std::cout << r << std::endl;
                    // still checks for EPSILON as theta could be as small as wanted.
                    if (r > EPSILON)
                    {
                        fx += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_x_center-x))/r;
                        fy += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_y_center-y))/r;
                    }
                    else
                    {
                        fx += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_x_center-x))/EPSILON;
                        fy += (G*m*tree[absOff + depthOff + i].m_av_mass*(tree[absOff + depthOff + i].m_y_center-y))/EPSILON;
                    }
                }
                else
                {
                    // recursive call
                    

                    absOff += nnode;
                    nnode*=4;
                    depthOff = (depthOff + i)*4;

                    compute_branch(absOff, depthOff, nnode, tree, x, y, m, fx, fy);
                }
            }
        }
    }
}
