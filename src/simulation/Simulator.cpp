#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args), numberofCells(0),
        logStep(0), saveStep(0), vlistStep(0), boxStep(0), box(0, 0, 0), 
        sb(params.output_file, params.surface_file, params.traj_file), traj(params.traj_file)
{
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
    drawBox = params.draw_box;
    pbc = params.pbc;
    gamma = params.k;
    Rc = params.r_cut;
    verlet_r = params.verlet_r;
    ttotal = params.ttime;
    nsteps = (int) ttotal / dt;
    setIntegrator(params.integrator_a);
    
    logStep = params.log_step;
    saveStep = params.save_step;
    boxStep = params.box_step; 
    vlistStep = params.vlist_step; 
    
    box.setX(params.bsx);
    box.setY(params.bsy);
    box.setZ(params.bsz);
    box.setDx(params.bsdx);
    box.setDy(params.bsdy);
    box.setDz(params.bsdz);
    box.setXend(0.5 * params.bsx);
    box.setYend(0.5 * params.bsy);
    box.setZend(0.5 * params.bsz);
    
    try
    {
        diagnoseParams();
    }
    catch(NotImplementedException& e)
    {
        cout <<  e.what() << endl;
        exit(1);
    }
    catch(DataException& e)
    {
        cout << e.what() << endl;
        exit(1);
    }
}

Simulator::Simulator(const Simulator& orig) : box(orig.box), sb(orig.sb), traj(orig.traj) {}

Simulator::~Simulator() {}

void Simulator::diagnoseParams()
{
    if (d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                "Single point representation is not implemented yet. "
                "Simulator is about to terminate !");
    if (d > 9)
        throw DataException("DataException:\n"
                "Depth of a triangulation to large ! "
                "For machine's safety Simulator is going to terminate !");
    
    if (pbc)
        throw NotImplementedException("NotImplementedException:\n"
                "Periodic boundary conditions are  not implemented yet."
                "Simulator is about to terminate !");
}

void Simulator::addCell(const Cell& newCell)
{
    try
    {
        if (cells.size() == MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached."
                                    "New cell will not be added !");
        
        cells.push_back(newCell);
        numberofCells++;        
    } 
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
}

void Simulator::addCell()
{
    addCell(1.0);
}

void Simulator::addCell(double r0)
{
    
    try
    {
        if (cells.size() == MAX_CELLS) 
            throw MaxSizeException("Maximum number of cells reached."
                                    "New cell will not be added !");
        
        SimpleTriangulation sm(params.d);
        list<Triangle> tris = sm.triangulate(r0);
        Cell newCell(tris);
   
        newCell.setA(a);
        newCell.setDp(dp);
        newCell.setRc(Rc);
        newCell.setRCellBox(params.r_bc);
        newCell.setGamma(gamma);
        newCell.setVerletR(verlet_r);
        newCell.setCellId(numberofCells);
        newCell.setMass(params.mass);
        newCell.setVisc(params.visc);
        newCell.setInitR(r0);
        
        addCell(newCell);       
    } 
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
}

void Simulator::initCells(int N, double r0)
{
    initCells(N, r0, r0);
}

void Simulator::initCells(int N, double ra, double rb)
{
    
    double nx, ny, nz;
    bool flag = true;
    double rbc = params.r_bc;
    double rsum;
    double r0;
    while(numberofCells < N)
    {
        r0 = uniform(ra, rb);
        flag = true;
        nx = uniform(-box.getX() + r0 + rbc, box.getX() - r0 - rbc);
        ny = uniform(-box.getY() + r0 + rbc, box.getY() - r0 - rbc);
        nz = uniform(-box.getZ() + r0 + rbc, box.getZ() - r0 - rbc);
        Vector3D shift(nx, ny, nz);
        for (int i = 0; i < numberofCells; i++)
        {
            Vector3D tmpcm = cells[i].getCm();
            Vector3D delta = shift - tmpcm;
            rsum = cells[i].getInitR() + r0 + Rc + EPSILON;
            if ( delta.length() < rsum)
            {
               flag = false; 
            }
        }
        if (flag)
        {
            addCell(r0);
            moveCell(shift, numberofCells-1);
        }
    }
}

void Simulator::simulate()
{
    simulate(nsteps);       
}

void Simulator::simulate(int steps)
{   
    rebuildVerletLists();
    sb.saveRenderScript(cells, box, drawBox);
    sb.saveSurfaceScript(cells);
    
    traj.open();
    traj.save(cells, getTotalVertices());

    for (int i = 0; i < steps; i++)
    {   
        if ( (i+1) % boxStep == 0)
        {
            box.resize();
        }
        
        integrate();
        if ( (i+1) % vlistStep == 0)
        {
            rebuildVerletLists();
        }
        
        if ( (i+1) % saveStep == 0)
        {
            traj.save(cells, getTotalVertices());
        }
    }
    traj.close();    
}

void Simulator::rebuildVerletLists()
{
    for (int i = 0; i < numberofCells; i++)
    {
        cells[i].voidVerletLsit();
    }
    
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < numberofCells; j ++)
        {
           cells[i].builtVerletList(cells[j]);    
        }
    }
}

void Simulator::calcForces()
{
    // RESET FORCES
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].voidForces();
    }
    
    // CALCULATE INTRA-CELLULAR FORCES 
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces();
    }
    
    // CALCULATE FORCES BETWEEN CELLS AND BOX
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcForces(box);
    }    
    
    // CALCULATE INTER-CELLULAR FORCES
    for (int i = 0; i < numberofCells; i++) 
    {
        for (int j = 0; j < numberofCells; j++) 
        {
            //if (i != j)
            //{
                //cells[i].calcForces(cells[j]);
                cells[i].calcForcesVL(cells[j]);   
            //}
        }
    }
}

void Simulator::setIntegrator(void (Simulator::*functoall)())
{
    integrator = functoall;
}

void Simulator::setIntegrator(char* token)
{
    if (STRCMP (token, "fe"))
    {
        this->setIntegrator(&Simulator::integrateEuler);
    }
    else if (STRCMP (token, "hm"))
    {
        this->setIntegrator(&Simulator::heunMethod);
    }
    else if (STRCMP (token, "rk"))
    {
        this->setIntegrator(&Simulator::midpointRungeKutta);
    }
    else if (STRCMP (token, "vv"))
    {
        this->setIntegrator(&Simulator::integrateVv);
    }
    else
    {
        this->setIntegrator(&Simulator::integrateEuler);
    }
}

void Simulator::integrate()
{
    (*this.*integrator)();
}

void Simulator::integrateEuler()
{
    calcForces();
    double m;
    double f;
    double mf;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            f = cells[i].vertices[j].getVisc();
            mf = m*f;
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / mf;
        }
    
    }
}

void Simulator::heunMethod()
{
    calcForces();
    double m;
    double f;
    double mf;
    
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
            cells[i].vertices[j].tmp_force = cells[i].vertices[j].force;
        }
    }
    
    //move the whole time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            f = cells[i].vertices[j].getVisc();
            mf = m*f;
            //cout << "m= "<< m << " f= "<< f << " mf= " << mf << endl;
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / mf;
        }
    }
    calcForces();
    
    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            f = cells[i].vertices[j].getVisc();
            mf = m*f;
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + 0.5 * dt * ( cells[i].vertices[j].tmp_force + cells[i].vertices[j].force) / mf;
        }
    }  
}

void Simulator::midpointRungeKutta()
{
    calcForces();
    double m;
    double f;
    double mf;
    
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
        }
    }
    
    //move half time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            f = cells[i].vertices[j].getVisc();
            mf = m*f;
            cells[i].vertices[j].xyz += 0.5 * dt * cells[i].vertices[j].force / mf;
        }
    }
    calcForces();
    
    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            f = cells[i].vertices[j].getVisc();
            mf = m*f;
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + dt * cells[i].vertices[j].force / mf;
        }
    }    
}

void Simulator::integrateVv()
{
    calcForces();
    double m;
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
        }
    }
    
    calcForces();
    for (int i = 0; i < numberofCells; i++) {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; // v(t+1) = v(t+1)_1 + 0.5 * dt * a(t+1)
        }
    }
}

void Simulator::moveCell(const Vector3D& v3d, int cellid)
{
    cells[cellid].addXYZ(v3d);
}

void Simulator::addCellVel(const Vector3D& v3d, int cellid)
{
    cells[cellid].addVelocity(v3d);
}

void Simulator::setBoxSize(const double bs)
{
    box.setX(bs);
    box.setY(bs);
    box.setZ(bs);
}

int Simulator::getTotalVertices()
{
    int totalnumber = 0;
    for (int i = 0; i < numberofCells; i++)
    {
        totalnumber += cells[i].numberofVertices();
    }
    return totalnumber;
}