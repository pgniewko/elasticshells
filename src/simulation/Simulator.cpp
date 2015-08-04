#include "Simulator.h"

utils::Logger Simulator::simulator_logs("simulator");

Simulator::Simulator(const arguments& args) : number_of_cells(0), box(0, 0, 0),
    sb(args.render_file, args.surface_file, args.traj_file, args.stress_file),
    traj(args.traj_file), log_sim(args.output_file, args.ob_config_file)
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

    params.log_step = args.log_step;
    params.save_step = args.save_step;
    params.box_step = args.box_step;
    params.vlist_step = args.vlist_step;
    params.d = args.d;
    params.nbhandler = args.nb_flag;
    params.E_cell = args.E_cell;
    params.nu = args.nu;
    params.th = args.thickness;
    params.dt = args.dt;
    params.dp = args.dp;
    params.ddp = args.ddp;
    params.visc = args.visc;
    params.mass = args.mass;
    params.ttime = args.ttime;
    params.r_vertex = args.r_vertex;
    params.verlet_r = args.verlet_r;
    params.growth_rate = args.growth_rate;
    params.vc = args.vc;
    params.bud_d = args.bud_d;
    params.div_ratio = args.div_ratio;
    params.draw_box = args.draw_box;
    params.scale = args.scale_flag;
    params.nsteps = args.nsteps ? args.nsteps : (int)params.ttime / params.dt;
    params.platotype = args.platotype;
    setIntegrator(args.integrator_a);
    setTriangulator(args.tritype);
    box.setX(args.bsx);
    box.setY(args.bsy);
    box.setZ(args.bsz);
    box.setDx(args.bsdx);
    box.setDy(args.bsdy);
    box.setDz(args.bsdz);
    box.setXstart(args.bsx);
    box.setYstart(args.bsy);
    box.setZstart(args.bsz);
    box.setXend(args.bsxe);
    box.setYend(args.bsye);
    box.setZend(args.bsze);
    box.setPbc(args.pbc);
    box.setEwall(args.E_wall);
    box.setNu(args.nu);
    domains.setupDomainsList(getMaxLengthScale(), box);
    OsmoticForce::setVolumeFlag(args.osmotic_flag);
    OsmoticForce::setEpsilon(args.eps);
    logParams();
}

Simulator::Simulator(const Simulator& orig) : number_of_cells(orig.number_of_cells),
    params(orig.params),
    box(orig.box), sb(orig.sb), traj(orig.traj),
    log_sim(orig.log_sim)
{
    exit(EXIT_FAILURE);
    // throw an exception - disallowed behavior
}

Simulator::~Simulator() {}

void Simulator::diagnoseParams(arguments args)
{
    if (args.d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                                      "Single point representation is not implemented yet. "
                                      "Simulator is about to terminate !");

    if (args.d > 7)
        throw DataException("DataException:\n"
                            "Depth of a triangulation too large ! "
                            "For machine's safety Simulator is about to terminate !");

    if (args.d < 0)
        throw DataException("Depth of a triangulation cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.E_wall < 0)
        throw DataException("Effective cell-box Young's modulus cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.E_cell < 0)
        throw DataException("Effective cell-cell Young's modulus cannot be negative\n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.dt <= 0)
        throw DataException("Time step must be positive number ! \n!"
                            "Simulation will terminate with exit(1)!\n");

//    if (args.k <= 0)
//        throw DataException("Spring constant for bonded vertices must be positive! \n!"
//                            "Simulation will terminate with exit(1)!\n");

    if (args.r_vertex <= 0)
        throw DataException("Vertex radius must be larger than 0! \n!"
                            "Simulation will terminate with exit(1)!\n");

    if (args.growth_rate < 0)
        throw DataException("Growth rate must be positive or 0! \n!"
                            "Simulation will terminate with exit(1)!\n");
}

void Simulator::logParams()
{
    simulator_logs << utils::LogLevel::INFO  << "SIM_STEPS=" << params.nsteps << "\n";
    simulator_logs << utils::LogLevel::INFO  << "BOX_STEP="  << params.box_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "SAVE_STEP=" << params.save_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "LOG_STEP="  << params.log_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "VERLET_STEP="  << params.vlist_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "TRIANGULATOR="  << triangulator << "\n";
    simulator_logs << utils::LogLevel::FINE  << "TIME STEP(DT)="  << params.dt << " [s]\n";
    simulator_logs << utils::LogLevel::FINE  << "DEPTH="  << params.d << "\n";
    simulator_logs << utils::LogLevel::FINE  << "DP="  << params.dp << " [bar]\n";
    simulator_logs << utils::LogLevel::FINE  << "DDP="  << params.ddp << " [bar]\n";
    simulator_logs << utils::LogLevel::FINE  << "E CELL="  << 0.1 * params.E_cell << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "E BOX="  << 0.1 * box.getE() << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "SURFACE_MODULUS="  << (params.E_cell * params.th) << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (CELL)="  << params.nu << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (BOX)="  << box.getNu() << "\n";
    simulator_logs << utils::LogLevel::FINE  << "R_VERTEX="  << params.r_vertex << " [micron]\n";
    simulator_logs << utils::LogLevel::FINE  << "GROWTH_RATE="  << params.growth_rate << "\n";
    simulator_logs << utils::LogLevel::FINE  << "VOLUME_CR="  << params.vc << "\n";
    simulator_logs << utils::LogLevel::FINE  << "BUD_SCAR_D="  << params.bud_d << "\n";
    simulator_logs << utils::LogLevel::FINE  << "CELL_DIV_RATIO="  << params.div_ratio << "\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.PBC=" << (box.pbc ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.BOX_DRAW=" << (params.draw_box ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_FLAG=" << (OsmoticForce::getFlag() ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_EPS=" << OsmoticForce::getEpsilon() << "\n";
    simulator_logs << utils::LogLevel::FINER << "MAX_SCALE=" << domains.getMaxScale() << "\n";
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

void Simulator::initCells(int N, double r0)
{
    initCells(N, r0, r0);
}

void Simulator::initCells(int N, double ra, double rb)
{
    if (ra > rb)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Illegal arguments: ra > rb. Simulator will set: rb = ra \n";
        rb = ra;
    }

    double nx, ny, nz;
    bool flag = true;
    double rc = 2.0 * params.r_vertex;
    double rsum;
    double r0;

    while (number_of_cells < N)
    {
        r0 = uniform(ra, rb);
        flag = true;
        nx = uniform(-box.getX() + r0 + rc, box.getX() - r0 - rc);
        ny = uniform(-box.getY() + r0 + rc, box.getY() - r0 - rc);
        nz = uniform(-box.getZ() + r0 + rc, box.getZ() - r0 - rc);
        
        if (number_of_cells == 0)
        {
            nx = 0;
            ny = 0;
            nz = 0;
        }

        Vector3D shift(nx, ny, nz);
         
        for (int i = 0; i < number_of_cells; i++)
        {
            Vector3D tmpcm = cells[i].getCm();
            Vector3D delta = shift - tmpcm;
            rsum = cells[i].getInitR() + r0 + rc + EPSILON;

            if ( delta.length() < rsum)
            {
                flag = false;
            }
        }

        if (flag)
        {
            addCell(r0);
            shiftCell(shift, number_of_cells - 1);
        }
    }
}

void Simulator::addCell(const Cell& newCell)
{

    try
    {
        if (cells.size() == MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached."
                                   "New cell will not be added !");
;
        cells.push_back(newCell);
        number_of_cells = cells.size();
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::WARNING << e.what() << "\n";
        return;
    }
}

void Simulator::addCell(double r0)
{
    try
    {
        if (cells.size() >= MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached.\n"
                                   "New cell will not be added !\n"
                                   "Program is going to TERMINANTE!");

        std::list<Triangle> tris;

        if (STRCMP (triangulator, "plato"))
        {
            PlatonicTriangulatoin tio(params.d, params.platotype);
            tris = tio.triangulate(r0);
        }
        else if (STRCMP (triangulator, "simple"))
        {
            SimpleTriangulation sm(params.d);
            tris = sm.triangulate(r0);
        }

        Cell newCell(tris);
        newCell.setEcc(params.E_cell);
        newCell.setNu(params.nu);
        newCell.setSpringConst(params.E_cell * params.th);
        newCell.setDp(params.dp, params.ddp);
        newCell.setVertexR(params.r_vertex);
        newCell.setVerletR(params.verlet_r);
        newCell.setCellId(number_of_cells);
        newCell.setMass(params.mass);
        newCell.setVisc(params.visc);
        newCell.setInitR(r0);
        newCell.setGrowthRate(params.growth_rate);
        newCell.setBuddingVolume(params.vc);
        newCell.setBudDiameter(params.bud_d);
        newCell.setDivisionRatio(params.div_ratio);
        addCell(newCell);
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}

void Simulator::simulate()
{
    simulate(params.nsteps);
}

void Simulator::simulate(int steps)
{
    if (params.nbhandler == 1)
    {
        rebuildVerletLists();
    }
    else if (params.nbhandler == 2)
    {
        rebuildDomainsList();
    }

    sb.saveRenderScript(cells, box, params.draw_box, 0.1);
    sb.saveSurfaceScript(cells);
    traj.open();
    traj.save(cells, getTotalVertices());
    log_sim.registerObservers();
    log_sim.open();
    log_sim.printHeader();
    log_sim.dumpState(box, cells);

    for (int i = 0; i <= steps; i++)
    {
        if (params.nbhandler == 1)
        {
            if ( (i + 1) % params.vlist_step == 0)
            {
                rebuildVerletLists();
            }
        }
        else if (params.nbhandler == 2)
        {
            rebuildDomainsList();
        }

        integrate();

        if ( (i + 1) % params.save_step == 0)
        {
            if (!params.scale)
            {
                traj.save(cells, getTotalVertices());
            }
            else
            {
                traj.save(cells, getTotalVertices(), box.getXstart() / box.getX(),
                          box.getYstart() / box.getY(), box.getZstart() / box.getZ());
            }
        }

        if ( (i + 1) % params.log_step == 0)
        {
            //log_sim.dumpState(box, cells, params.r_vertex, (i + 1), getTotalVertices(), params.nbhandler);
            log_sim.dumpState(box, cells);
        }

        if ( (i + 1) % params.box_step == 0)
        {
            box.resize();
            domains.setBoxDim(box);
        }

        for (int i = 0; i < number_of_cells; i++)
        {
            cells[i].cellCycle(params.dt);
        }

        if ( i % (steps / 10) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * i / steps << "% OF THE SIMULATION IS DONE" "\n";
        }
    }

    sb.saveStressScript(cells, box);
    traj.close();
    log_sim.close();
}

void Simulator::calcForces()
{
    // RESET FORCES
    for (int i = 0 ; i < number_of_cells; i++)
    {
        cells[i].voidForces();
    }

    // CALCULATE INTRA-CELLULAR FORCES
    for (int i = 0 ; i < number_of_cells; i++)
    {
        cells[i].calcBondedForces();
    }

    // CALCULATE INTER-CELLULAR FORCES
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < number_of_cells; j++)
        {
            if (params.nbhandler == 0)
            {
                cells[i].calcNbForcesON2(cells[j], box);
            }
            else if (params.nbhandler == 1)
            {
                cells[i].calcNbForcesVL(cells[j], box);
            }
            else if (params.nbhandler == 2)
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
        for (int i = 0 ; i < number_of_cells; i++)
        {
            cells[i].calcBoxForces(box);
        }
    }
}

void Simulator::rebuildVerletLists()
{
    for (int i = 0; i < number_of_cells; i++)
    {
        cells[i].voidVerletLsit();
    }

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < number_of_cells; j ++)
        {
            cells[i].builtVerletList(cells[j], box);
        }
    }
}

void Simulator::rebuildDomainsList()
{
    domains.voidDomains();

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            domains.assignVertex(cells[i].vertices[j], i);
        }
    }

    for (int i = 0; i < number_of_cells; i++)
    {
        cells[i].voidVerletLsit();
    }

    for (int i = 0; i < number_of_cells; i++)
    {
        cells[i].builtNbList(cells, domains, box);
    }
}

void Simulator::shiftCell(const Vector3D& v3d, int cellid)
{
    cells[cellid].addXYZ(v3d);
}

int Simulator::getTotalVertices()
{
    int totalnumber = 0;

    for (int i = 0; i < number_of_cells; i++)
    {
        totalnumber += cells[i].getNumberVertices();
    }

    return totalnumber;
}

double Simulator::getMaxLengthScale()
{
    double maxscale = 0.0;
    maxscale = std::max(maxscale, params.r_vertex);
    return maxscale;
}

/*
 * TRIANGULATION
 */

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
        simulator_logs << utils::LogLevel::FINE  << "DEFAULT FORWARD EULER IS USED\n";
    }
}

void Simulator::setTriangulator(char* token)
{
    if (STRCMP (token, "simple"))
    {
        triangulator = token;
    }
    else if (STRCMP (token, "plato"))
    {
        triangulator = token;
    }
    else
    {
        triangulator = (char*)& "simple";
        simulator_logs << utils::LogLevel::FINE  << "SIMPLE TRIANGULATION IS APPLIED\n";
    }
}

/*
 * INTEGRATORS
 */

void Simulator::integrate()
{
    (*this.*integrator)();
//    makeVertsOlder(); // function's temporary location
}

void Simulator::setIntegrator(void (Simulator::*functoall)())
{
    integrator = functoall;
}

void Simulator::integrateEuler()
{
    calcForces();
    double visc;
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
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
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
            cells[i].vertices[j].tmp_force = cells[i].vertices[j].force;
        }
    }

    //move the whole time-step and calculate  forces
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].force / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
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
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].tmp_xyz = cells[i].vertices[j].xyz;
        }
    }

    //move half time-step and calculate  forces
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].xyz += 0.5 * dt * cells[i].vertices[j].force / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
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
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].xyz += dt * cells[i].vertices[j].velocity; // x(t+1)_a = x(t) + v(t)*dt
            cells[i].vertices[j].xyz += 0.5 * dt * dt * cells[i].vertices[j].force / m; // x(t+1) = x(t+1)_a + 0.5*dt*dt* a(t)
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; //v(t+1)_1 = v(t) + 0.5 * dt * a(t)
        }
    }

    calcForces();

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            m = cells[i].vertices[j].getMass();
            cells[i].vertices[j].velocity += 0.5 * dt * cells[i].vertices[j].force / m; // v(t+1) = v(t+1)_1 + 0.5 * dt * a(t+1)
        }
    }
}