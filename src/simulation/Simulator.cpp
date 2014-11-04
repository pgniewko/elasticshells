#include "Simulator.h"

Simulator::Simulator(const arguments& args) : params(args), numberofCells(0),
        dt(0), a(0), dp(0), gamma(0), R0(0), Rc(0), ttotal(0), initcellmass(0),
        verlet_r(0), nsteps(0), d(0),
        logStep(1), saveStep(1), vlistStep(1), boxStep(1), 
        box(0, 0, 0), drawBox(false), nbhandler(0),
        sb(params.render_file, params.surface_file, params.traj_file), 
        traj(params.traj_file), 
        logsim(params.output_file),
        simulator_logs("simulator")
{
    
    try
    {
        diagnoseParams(args);
    }
    catch (NotImplementedException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    catch (DataException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    
    dt = params.dt;
    a = params.a;
    d = params.d;
    dp = params.dp;
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
    box.setXend(params.bsxe);
    box.setYend(params.bsye);
    box.setZend(params.bsze);
    drawBox = params.draw_box;
    box.setPbc(params.pbc);
    box.setEcw(params.ecw);
    nbhandler = params.nbFlag;

    double maxscale = getMaxScale();
    domains.setupDomainsList(maxscale, box);
    OsmoticForce::setVolumeFlag(params.osmFlag);
    
    logParams();
}

Simulator::Simulator(const Simulator& orig) : params(orig.params), numberofCells(orig.numberofCells),
        dt(orig.dt), a(orig.a), dp(orig.dp), gamma(orig.gamma), R0(orig.R0), 
        Rc(orig.Rc), ttotal(orig.ttotal), initcellmass(orig.initcellmass),
        verlet_r(orig.verlet_r), nsteps(orig.nsteps), d(orig.d),
        logStep(orig.logStep), saveStep(orig.saveStep), vlistStep(orig.vlistStep), boxStep(orig.boxStep), 
        box(orig.box), drawBox(false), nbhandler(orig.nbhandler),
        sb(orig.sb), traj(orig.traj), logsim(orig.logsim), simulator_logs(orig.simulator_logs)
{
    // exception disallowed behavior
}

Simulator::~Simulator() {}

void Simulator::diagnoseParams(arguments args)
{
    if (args.d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                                      "Single point representation is not implemented yet. "
                                      "Simulator is about to terminate !");

    if (args.d > 9)
        throw DataException("DataException:\n"
                            "Depth of a triangulation too large ! "
                            "For machine's safety Simulator is about to terminate !");
    
    if (args.d < 0)
        throw DataException("Depth of a triangulation cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");
}

void Simulator::logParams()
{
    
    simulator_logs << utils::LogLevel::INFO << "BOX_STEP="  << boxStep << "\n";
    simulator_logs << utils::LogLevel::INFO << "SAVE_STEP=" << saveStep << "\n";
    simulator_logs << utils::LogLevel::INFO << "LOG_STEP="  << logStep << "\n";
    simulator_logs << utils::LogLevel::INFO << "VERLET_STEP="  << vlistStep << "\n";
    
    simulator_logs << utils::LogLevel::FINE << "TIME STEP(DT)="  << dt << "\n";
    simulator_logs << utils::LogLevel::FINE << "DEPTH="  << d << "\n";
    simulator_logs << utils::LogLevel::FINE << "DP="  << dp << "\n";
    simulator_logs << utils::LogLevel::FINE << "GAMMA="  << gamma << "\n";
    simulator_logs << utils::LogLevel::FINE << "A="  << a << "\n";
    simulator_logs << utils::LogLevel::FINE << "E* CELL-BOX="  << box.ecw << "\n";
    simulator_logs << utils::LogLevel::FINE << "R:CELL_CELL="  << Rc << "\n";
    simulator_logs << utils::LogLevel::FINE << "R:CELL_BOX="  << params.r_bc << "\n";
    simulator_logs << utils::LogLevel::FINE << "BOX.PBC="<<(box.pbc ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINE << "BOX.BOX_DRAW=" << (drawBox ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINE << "OSMOTIC_FLAG=" << (params.osmFlag ? "true" : "false") << "\n";
    
    simulator_logs << utils::LogLevel::FINER << "BOX.X="  << box.getX() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.X="  << box.getY() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Z="  << box.getZ() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.DX=" << box.getDx() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.DX=" << box.getDy() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.DZ=" << box.getDz() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.XE=" << box.getXend() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.YE=" << box.getYend() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.ZE=" << box.getZend() << "\n";
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
        simulator_logs << utils::LogLevel::WARNING << e.what() << "\n";
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
        if (cells.size() >= MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached.\n"
                                   "New cell will not be added !\n"
                                   "Program is going to TERMINANTE!");

        SimpleTriangulation sm(params.d);
        std::list<Triangle> tris = sm.triangulate(r0);
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
        newCell.setNRT(params.dp);
        addCell(newCell);
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::WARNING << e.what() << "\n";
        exit(1);
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

    while (numberofCells < N)
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
            moveCell(shift, numberofCells - 1);
        }
    }
}

void Simulator::simulate()
{
    simulate(nsteps);
}

void Simulator::simulate(int steps)
{
    if (nbhandler == 1)
    {
        rebuildVerletLists();
    }
    else if (nbhandler == 2)
    {
        rebuildDomainsList();
    }
    
    sb.saveRenderScript(cells, box, drawBox);
    sb.saveSurfaceScript(cells);
    traj.open();
    traj.save(cells, getTotalVertices());

    logsim.open();
    logsim.dumpState(box, cells, 1, getTotalVertices(), nbhandler);
    
    
    for (int i = 0; i <= steps; i++)
    {
        if (nbhandler == 1)
        {
            if ( (i + 1) % vlistStep == 0)
            {
                rebuildVerletLists();
            }
        }
        else if (nbhandler == 2)
        {
            rebuildDomainsList();
        }
        
        integrate();

        if ( (i + 1) % saveStep == 0)
        {
            if (false)
            {
               traj.save(cells, getTotalVertices());
            }
            else
            {
                traj.save(cells, getTotalVertices(), params.bsx / box.getX(), params.bsy / box.getY(), params.bsz / box.getZ());
            }
        }
        
        if ( (i+1) % logStep == 0)
        {
            logsim.dumpState(box, cells, (i+1), getTotalVertices(), nbhandler);
        }
        
        if ( (i + 1) % boxStep == 0)
        {
            box.resize();
            domains.setBoxDim(box);
        }
        
        if ( i % (steps / 10) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * i / steps << "% OF THE SIMULATION IS DONE" "\n";
        }
    }

    
    traj.close();
    logsim.close();
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
        cells[i].calcBondedForces();
    }

    // CALCULATE INTER-CELLULAR FORCES
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < numberofCells; j++)
        {
            if (nbhandler == 0)
            {
                cells[i].calcNbForcesON2(cells[j], box);
            }
            else if (nbhandler == 1)
            {
                cells[i].calcNbForcesVL(cells[j], box);
            }
            else if (nbhandler == 2)
            {
                cells[i].calcNbForcesVL(cells[j], box);
            }
            else 
            {
                cells[i].calcNbForcesON2(cells[j], box);
            }
        }
    }    
    
    // CALCULATE FORCES BETWEEN CELLS AND BOX
    if (!box.pbc) 
    {
        for (int i = 0 ; i < numberofCells; i++)
        {
            cells[i].calcBoxForces(box);
        }
    }
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
            cells[i].builtVerletList(cells[j], box);
        }
    }
}

void Simulator::rebuildDomainsList()
{
    
    domains.voidDomains();
    
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
             domains.assignVertex(cells[i].vertices[j], i);
        }
       
    }
    
    for (int i = 0; i < numberofCells; i++)
    {
        cells[i].voidVerletLsit();
    }
    
    for (int i = 0; i < numberofCells; i++)
    {
        cells[i].builtNbList(cells, domains, box);
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
    double visc;

    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / visc;
        }
    }
}

void Simulator::heunMethod()
{
    calcForces();
    double visc;

    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
            cells[i].vertices[j].tmp_force = cells[i].vertices[j].force;
        }
    }

    //move the whole time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + 0.5 * dt * ( cells[i].vertices[j].tmp_force + cells[i].vertices[j].force) / visc;
        }
    }
}

void Simulator::midpointRungeKutta()
{
    calcForces();
    double visc;

    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
        }
    }

    //move half time-step and calculate  forces
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz += 0.5 * dt * cells[i].vertices[j].force / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz = cells[i].vertices[j].tmp_xyz + dt * cells[i].vertices[j].force / visc;
        }
    }
}

void Simulator::integrateVv()
{
    calcForces();
    double m;

    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
        }
    }

    calcForces();

    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
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

int Simulator::getNumberOfCells()
{
    return numberofCells;
}

int Simulator::getTotalVertices()
{
    int totalnumber = 0;

    for (int i = 0; i < numberofCells; i++)
    {
        totalnumber += cells[i].numberOfVerts();
    }

    return totalnumber;
}

double Simulator::getMaxScale()
{
    double maxscale = 0.0;
    maxscale = std::max(params.r_bc, Rc);
    return maxscale;
}