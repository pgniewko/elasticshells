#include "Simulator.h"

utils::Logger Simulator::simulator_logs("simulator");
ulong Simulator::FORCE_EVALUATION_COUTER(0);

Simulator::Simulator(const arguments& args) : number_of_cells(0), box(0, 0, 0),
    sb(args.render_file, args.surface_file, args.traj_file, args.stress_file),
    traj(args.traj_file, args.box_file), log_sim(args.output_file, args.ob_config_file),
    restarter(args.topology_file, args.lf_file)
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
    params.d = args.d;
    params.nbhandler = args.nb_flag;
    params.E_cell = args.E_cell;
    params.nu = args.nu;
    params.th = args.thickness;
    params.dt = args.dt;
    params.dp = args.dp;
    params.ddp = args.ddp;
    params.volume_scale = args.volume_scale;
    params.ttime = args.ttime;
    params.r_vertex = args.r_vertex;
    params.draw_box = args.draw_box;
    params.scale = args.scale_flag;
    params.dynamics = args.dynamics;
    params.const_volume = args.const_volume;
    params.nsteps = args.nsteps ? args.nsteps : (int)params.ttime / params.dt;
    params.platotype = args.platotype;
    setIntegrator(args.integrator_a);
    setTriangulator(args.tritype);
    box.setX(args.bsx);
    box.setY(args.bsy);
    box.setZ(args.bsz);
    box.setXmax(args.bsx);
    box.setYmax(args.bsy);
    box.setZmax(args.bsz);
    box.setXmin(args.bsxe);
    box.setYmin(args.bsye);
    box.setZmin(args.bsze);
    box.setPbc(args.pbc);
    box.setEwall(args.E_wall);
    box.setNu(args.nu);
    box.setDefaultSchedule(params.nsteps, args.box_step, args.bsdx, args.bsdy, args.bsdz, 0.0, 0.0, 0.0);
    box.configureScheduler(args.sch_config_file);

    domains.setupDomainsList(getLengthScale(), box);
    OsmoticForce::setVolumeFlag(args.osmotic_flag);
    OsmoticForce::setEpsilon(args.eps);
    Cell::no_bending = args.nobending;
    logParams();

}

Simulator::~Simulator() {}

void Simulator::diagnoseParams(arguments args)
{
    if (args.d == 0)
        throw NotImplementedException("NotImplementedException:\n"
                                      "Single point representation is not implemented yet. "
                                      "Simulator is about to terminate !");

    if (args.d > 8)
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

    if (args.r_vertex <= 0)
        throw DataException("Vertex radius must be larger than 0! \n!"
                            "Simulation will terminate with exit(1)!\n");

}

void Simulator::logParams()
{
    simulator_logs << utils::LogLevel::INFO  << "SIM_STEPS=" << params.nsteps << "\n";
    simulator_logs << utils::LogLevel::INFO  << "SAVE_STEP=" << params.save_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "LOG_STEP="  << params.log_step << "\n";
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
    simulator_logs << utils::LogLevel::FINE  << "BOX.PBC=" << (box.pbc ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.BOX_DRAW=" << (params.draw_box ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_FLAG=" << (OsmoticForce::getFlag() ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_EPS=" << OsmoticForce::getEpsilon() << "\n";
    simulator_logs << utils::LogLevel::FINER << "MAX_SCALE=" << domains.getMaxScale() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.X="  << box.getX() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Y="  << box.getY() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Z="  << box.getZ() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.XE=" << box.getXmin() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.YE=" << box.getYmin() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.ZE=" << box.getZmin() << "\n";
}

void Simulator::initCells(int N, double r0)
{
    initCells(N, r0, r0);
}

void Simulator::initCells(int N, double ra, double rb)
{
    initCells(N, ra, rb, (char*)&"ms_kot", false);
}

void Simulator::initCells(int N, double ra, double rb, char* model_t, bool restart_flag)
{
    if (ra > rb)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Illegal arguments: ra > rb. Simulator will set: rb = ra \n";
        rb = ra;
    }

    if (restart_flag)
    {
        simulator_logs << utils::LogLevel::INFO  << "Cells are initialized from the topology file\n";
    }

    simulator_logs << utils::LogLevel::INFO  << "CELL MODEL: " << model_t << "\n";

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
            rsum = cells[i].getInitR() + r0 + rc + constants::epsilon;

            if ( delta.length() < rsum)
            {
                flag = false;
            }
        }

        if (flag)
        {
            addCell(r0, model_t);
            shiftCell(shift, number_of_cells - 1);
        }
    }

    restarter.saveTopologyFile(cells, model_t);
    set_min_force();
}

void Simulator::pushCell(const Cell& newCell)
{
    try
    {
        if (cells.size() == MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached."
                                   "New cell will not be added !");

        cells.push_back(newCell);
        number_of_cells = cells.size();
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::WARNING << e.what() << "\n";
        return;
    }
}

void Simulator::addCell(double r0, char* model_t)
{
    try
    {
        if (cells.size() >= MAX_CELLS)
            throw MaxSizeException("Maximum number of cells reached.\n"
                                   "New cell will not be added !\n"
                                   "Program is going to TERMINANTE!");

        std::list<Triangle> tris;

        if ( !triangulator.compare("plato") )
        {
            PlatonicTriangulatoin tio(params.d, params.platotype);
            tris = tio.triangulate(r0);
        }
        else if ( !triangulator.compare("simple") )
        {
            SimpleTriangulation sm(params.d);
            tris = sm.triangulate(r0);
        }

        Cell newCell(tris);

        newCell.setEcc(params.E_cell);
        newCell.setNu(params.nu);
        newCell.setSpringConst(params.E_cell, params.th, params.nu, model_t);
        newCell.setBSprings(params.E_cell, params.th, params.nu);
        newCell.setDp(params.dp, params.ddp);
        newCell.setConstantVolume(params.volume_scale);
        newCell.setVertexR(params.r_vertex);
        newCell.setCellId(number_of_cells);
        newCell.setInitR(r0);

        pushCell(newCell);
    }
    catch (MaxSizeException& e)
    {
        simulator_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}

void Simulator::restart()
{
    restarter.readTopologyFile(cells);
}


void Simulator::simulate()
{
    simulate(params.nsteps);
}

void Simulator::simulate(int steps)
{
    updateCells();

//    if (params.nbhandler == 1)
//    {
//        rebuildVerletLists();
//    }
    if (params.nbhandler == 2)
    {
        rebuildDomainsList();
    }

    sb.saveRenderScript(cells, box, params.draw_box, 0.1);
    sb.saveSurfaceScript(cells);
    traj.open();
    traj.save_traj(cells, getTotalVertices());
    log_sim.registerObservers();
    log_sim.open();
    log_sim.printHeader();
    log_sim.dumpState(box, cells);

    bool resized = false;

    //calcForces();
    for (int i = 0; i < steps; i++)
    {

        if ( i % (steps / std::min(steps, 10) ) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * i / steps << "% OF THE SIMULATION IS DONE" "\n";
        }

        do
        {
            do
            {
                integrate();
            }
            while ( check_min_force() );
        }
        while ( check_const_volume() );


        if ( (i + 1) % params.save_step == 0)
        {
            if (!params.scale)
            {
                traj.save_traj(cells, getTotalVertices());
            }
            else
            {
                traj.save_traj(cells, getTotalVertices(), box.getXmax() / box.getX(),
                               box.getYmax() / box.getY(), box.getZmax() / box.getZ());
            }

            traj.save_box(box, (i + 1) * params.dt);
        }

        if ( (i + 1) % params.log_step == 0)
        {
            log_sim.dumpState(box, cells);
            restarter.saveLastFrame(cells);
        }

        resized = box.resize();

        if (resized)
        {
            domains.setBoxDim(box);
        }
    }

    log_sim.dumpState(box, cells); // TODO: fix that. the forces are not updated etc. That's causing weird results, probably there is not force relaxation before dump
    restarter.saveLastFrame(cells);
    sb.saveStrainScript(cells, box);
    traj.close();
    log_sim.close();

    simulator_logs << utils::LogLevel::FINEST << "Forces have been evaluated: " << FORCE_EVALUATION_COUTER << " times.\n";
    simulator_logs << utils::LogLevel::FINEST << "Energy  has been evaluated: " << Energy::ENERGY_EVALUATION_COUNTER << " times.\n";
}

void Simulator::calcForces()
{
    FORCE_EVALUATION_COUTER++;
    #pragma omp parallel
    {
        // CALC CENTER OF MASS
        #pragma omp for
        for (uint i = 0; i < cells.size(); i++)
        {
            cells[i].update();
        }

        // RESET FORCES
        #pragma omp for

        for (int i = 0 ; i < number_of_cells; i++)
        {
            cells[i].voidForces();
        }

        // CALCULATE INTRA-CELLULAR FORCES
        #pragma omp for schedule(guided)

        for (int i = 0 ; i < number_of_cells; i++)
        {
            cells[i].calcBondedForces();
        }

        // CALCULATE INTER-CELLULAR FORCES
        if (params.nbhandler == 0)
        {
            #pragma omp for schedule(guided)

            for (int i = 0; i < number_of_cells; i++)
            {
                for (int j = 0; j < number_of_cells; j++)
                {
                    cells[i].calcNbForcesON2(cells[j], box);
                }
            }
        }
        else if (params.nbhandler == 2)
        {
            domains.calcNbForces(cells, box);
        }

        // CALCULATE FORCES BETWEEN CELLS AND BOX
        if (!box.pbc)
        {
            #pragma omp for schedule(guided)

            for (int i = 0 ; i < number_of_cells; i++)
            {
                cells[i].calcBoxForces(box);
            }
        }
    }
}
void Simulator::update_neighbors_list()
{
    if (params.nbhandler == 2)
    {
        rebuildDomainsList();
    }
}

void Simulator::rebuildDomainsList()
{
    domains.voidDomains();

    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            domains.assignVertex(&cells[i].vertices[j]);
        }
    }
}

void Simulator::shiftCell(const Vector3D& v3d, int cellid)
{
    cells[cellid].addXYZ(v3d);
    cells[cellid].update();
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

double Simulator::getLengthScale()
{
    double maxscale = 0.0;

    if (params.nbhandler == 2)
    {
        maxscale = std::max(maxscale, params.r_vertex);
    }

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
    else if (STRCMP (token, "cp"))
    {
        this->setIntegrator(&Simulator::gear_cp);
    }
    else if (STRCMP (token, "cg"))
    {
        this->setIntegrator(&Simulator::cg);
    }
    else
    {
        this->setIntegrator(&Simulator::integrateEuler);
        simulator_logs << utils::LogLevel::FINE  << "DEFAULT FORWARD EULER IS USED\n";
    }

    if (integrator == NULL)
    {
        simulator_logs << utils::LogLevel::CRITICAL << "integrator == NULL";
        exit(EXIT_FAILURE);
    }
}

void Simulator::setTriangulator(char* token)
{
    if (STRCMP (token, "simple"))
    {
        triangulator = std::string(token);
    }
    else if (STRCMP (token, "plato"))
    {
        triangulator = std::string(token);
    }
    else
    {
        triangulator = std::string("simple");
        simulator_logs << utils::LogLevel::FINE  << "SIMPLE TRIANGULATION IS APPLIED\n";
    }
}

void Simulator::updateCells()
{
    for (uint i = 0; i < cells.size(); i++)
    {
        cells[i].update();
    }
}

void Simulator::set_min_force()
{
    double average_area = cells[0].calcSurfaceArea();
    average_area /= cells[0].getNumberTriangles();

    double max_turgor = 0.0;

    for (int i = 0; i < number_of_cells; i++)
    {
        max_turgor = std::max(max_turgor, cells[i].getTurgor());
    }

    MIN_FORCE_SQ = FORCE_FRAC * max_turgor * average_area;
    MIN_FORCE_SQ = MIN_FORCE_SQ * MIN_FORCE_SQ;
}

bool Simulator::check_min_force()
{
    if (params.dynamics)
    {
        return false;
    }

    if (integrator == &Simulator::cg)
    {
        return false;
    }

    
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            if (cells[i].vertices[j].f_c.length_sq() > MIN_FORCE_SQ)
            {
                return true;
            }
        }
    }
    
    return false;
}

bool Simulator::check_const_volume()
{
    if (params.const_volume)
    {
        double step;
        double eps = 0.001;
        bool flag = false;

        for (int i = 0; i < number_of_cells; i++)
        {
            step = cells[i].checkVolumeCondition(eps);

            if ( fabs(step) > eps )
            {
                flag = true;
                cells[i].ajustTurgor(step);
            }
        }

        return flag;
    }
    else
    {
        return false;
    }
}

/*
 * INTEGRATORS
 *
 * Viscosity of each vertex is assumed to be 1.0 !
 *
 */
void Simulator::integrate()
{
    update_neighbors_list();
    (*this.*integrator)();
    updateCells();
}

void Simulator::setIntegrator(void (Simulator::*functoall)())
{
    integrator = functoall;
}

void Simulator::integrateEuler()
{
    calcForces();
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c += dt * cells[i].vertices[j].f_c;
        }
    }
}

void Simulator::heunMethod()
{
    calcForces();
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_p = cells[i].vertices[j].r_c;
            cells[i].vertices[j].f_p = cells[i].vertices[j].f_c;
        }
    }

    //move the whole time-step and calculate  forces
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c += dt * cells[i].vertices[j].f_c;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c = cells[i].vertices[j].r_p + 0.5 * dt * ( cells[i].vertices[j].f_p + cells[i].vertices[j].f_c);
        }
    }
}

void Simulator::midpointRungeKutta()
{
    double dt = params.dt;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_p = cells[i].vertices[j].r_c;
        }
    }

    calcForces();

    //move half time-step and calculate  forces
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c += 0.5 * dt * cells[i].vertices[j].f_c;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c = cells[i].vertices[j].r_p + dt * cells[i].vertices[j].f_c;
        }
    }
}

void Simulator::gear_cp()
{
    double dt = params.dt;
    double C1, C2;

    C1 = dt;
    C2 = dt * dt / 2.0;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_p = cells[i].vertices[j].r_c + C1 * cells[i].vertices[j].v_c + C2 * cells[i].vertices[j].a_c;
            cells[i].vertices[j].v_p = cells[i].vertices[j].v_c + C1 * cells[i].vertices[j].a_c;
            cells[i].vertices[j].a_p = cells[i].vertices[j].a_c;
        }
    }

    calcForces();

    double gear0 = 5.0 / 12.0;
    double gear2 = 1.0 / 2.0;

    double CR = gear0 * C1;
    double CA = gear2 * C1 / C2;

    Vector3D corr_v;

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            corr_v = cells[i].vertices[j].f_c - cells[i].vertices[j].v_p; // viscosity = 1.0

            cells[i].vertices[j].r_c = cells[i].vertices[j].r_p + CR * corr_v;
            cells[i].vertices[j].v_c = cells[i].vertices[j].f_c;
            cells[i].vertices[j].a_c = cells[i].vertices[j].a_p + CA * corr_v;
        }
    }
}

/*
 * **************************************
 * CONJUGATE GRADIENTS CODE STARTS HERE *
 * **************************************
 */
void Simulator::cg()
{
    int n = 3 * getTotalVertices();
    
    double* p = darray(n);
    
    int counter = 0;
    
    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            p[3 * counter + 0] = cells[i].vertices[j].r_c.x;
            p[3 * counter + 1] = cells[i].vertices[j].r_c.y;
            p[3 * counter + 2] = cells[i].vertices[j].r_c.z;
            counter++;
        }
    }

    double ftol = 1e-10;
    int iter;
    double fret;

    frprmn(p, n, ftol, &iter, &fret);
    counter = 0;
    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c.x = p[3 * counter + 0];
            cells[i].vertices[j].r_c.y = p[3 * counter + 1];
            cells[i].vertices[j].r_c.z = p[3 * counter + 2];
            counter++;
        }
    }
}

static double maxarg1, maxarg2;

double Simulator::func(double _p[])
{
    int counter = 0;

    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c.x = _p[3 * counter + 0];
            cells[i].vertices[j].r_c.y = _p[3 * counter + 1];
            cells[i].vertices[j].r_c.z = _p[3 * counter + 2];
            counter++;
        }
    }
    
    return Energy::calcTotalEnergy(cells, box, domains);
}

void Simulator::dfunc(double _p[], double _xi[])
{
    int counter = 0;

    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].r_c.x = _p[3 * counter + 0];
            cells[i].vertices[j].r_c.y = _p[3 * counter + 1];
            cells[i].vertices[j].r_c.z = _p[3 * counter + 2];
            counter++;
        }
    }

    calcForces();

    counter = 0;

    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            _xi[3 * counter + 0] = -cells[i].vertices[j].f_c.x;
            _xi[3 * counter + 1] = -cells[i].vertices[j].f_c.y;
            _xi[3 * counter + 2] = -cells[i].vertices[j].f_c.z;
            counter++;
        }
    }
}


#define ITMAX 100000
#define EPS 1.0e-10
#define FREEALL free_darray(xi); free_darray(h); free_darray(g);
void Simulator::frprmn(double p[], int n, double ftol, int* iter, double* fret)
{
    int j, its;
    double gg, gam, fp, dgg;
    double* g, *h, *xi;

    g = darray(n);
    h = darray(n);
    xi = darray(n);
    fp = func(p);
    dfunc(p, xi);

    for (j = 0; j < n; j++)
    {
        g[j] = -xi[j];
        xi[j] = h[j] = g[j];
    }

    for (its = 1; its <= ITMAX; its++)
    {
        *iter = its;
        linmin(p, xi, n, fret);
        //dlinmin(p, xi, n, fret);
       
        if (2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + EPS))
        {
            FREEALL
            return;
        }

        fp = *fret;
        dfunc(p, xi);
        dgg = gg = 0.0;

        for (j = 0; j < n; j++)
        {
            gg += g[j] * g[j];
            dgg += (xi[j] + g[j]) * xi[j];
        }

        if (gg == 0.0)
        {
            FREEALL
            return;
        }

        gam = dgg / gg;

        for (j = 0; j < n; j++)
        {
            g[j] = -xi[j];
            xi[j] = h[j] = g[j] + gam * h[j];
        }
    }

    nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL


#define TOL 1.0e-4
static int ncom;
static double* pcom;
static double* xicom;

void Simulator::linmin(double p[], double xi[], int n, double* fret)
{
    int j;
    double xx, xmin, fx, fb, fa, bx, ax;

    ncom = n;
    pcom = darray(n);
    xicom = darray(n);

    for (j = 0; j < n; j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }

    ax = 0.0;
    xx = 1.0;

    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb);
    *fret = brent(ax, xx, bx, TOL, &xmin);

    for (j = 0; j < n; j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }

    free_darray(xicom);
    free_darray(pcom);
}
#undef TOL

#define TOL 1.0e-4
void Simulator::dlinmin(double p[], double xi[], int n, double* fret)
{
    int j;
    double xx, xmin, fx, fb, fa, bx, ax;

    ncom = n;
    pcom = darray(n);
    xicom = darray(n);

    for (j = 0; j < n; j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }

    ax = 0.0;
    xx = 1.0;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb);
    *fret = dbrent(ax, xx, bx, TOL, &xmin);

    for (j = 0; j < n; j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }

    free_darray(xicom);
    free_darray(pcom);
}
#undef TOL

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void Simulator::mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc)
{
    double ulim, u, r, q, fu, dum;

    *fa = f1dim(*ax);
    *fb = f1dim(*bx);

    if (*fb > *fa)
    {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    }

    *cx = (*bx) + GOLD * (*bx - *ax);
    *fc = f1dim(*cx);

    while (*fb > *fc)
    {
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
            (2.0 * SIGN2(FMAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);

        if ((*bx - u) * (u - *cx) > 0.0)
        {
            fu = f1dim(u);

            if (fu < *fc)
            {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx = u;
                *fc = fu;
                return;
            }

            u = (*cx) + GOLD * (*cx - *bx);
            fu = f1dim(u);
        }
        else if ((*cx - u) * (u - ulim) > 0.0)
        {
            fu = f1dim(u);

            if (fu < *fc)
            {
                SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
                SHFT(*fb, *fc, fu, f1dim(u))
            }
        }
        else if ((u - ulim) * (ulim - *cx) >= 0.0)
        {
            u = ulim;
            fu = f1dim(u);
        }
        else
        {
            u = (*cx) + GOLD * (*cx - *bx);
            fu = f1dim(u);
        }

        SHFT(*ax, *bx, *cx, u)
        SHFT(*fa, *fb, *fc, fu)
    }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT


#define ITMAX 10000
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d);
double Simulator::brent(double ax, double bx, double cx, double tol, double* xmin)
{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = f1dim(x);

    for (iter = 1; iter <= ITMAX; iter++)
    {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);

        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
        {
            *xmin = x;
            return fx;
        }

        if (fabs(e) > tol1)
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);

            if (q > 0.0)
            {
                p = -p;
            }

            q = fabs(q);
            etemp = e;
            e = d;

            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
            {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            }
            else
            {
                d = p / q;
                u = x + d;

                if (u - a < tol2 || b - u < tol2)
                {
                    d = SIGN2(tol1, xm - x);
                }
            }
        }
        else
        {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }

        u = (fabs(d) >= tol1 ? x + d : x + SIGN2(tol1, d));
        fu = f1dim(u);

        if (fu <= fx)
        {
            if (u >= x)
            {
                a = x;
            }
            else
            {
                b = x;
            }

            SHFT(v, w, x, u)
            SHFT(fv, fw, fx, fu)
        }
        else
        {
            if (u < x)
            {
                a = u;
            }
            else
            {
                b = u;
            }

            if (fu <= fw || w == x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }

    nrerror("Too many iterations in brent");
    *xmin = x;
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define ITMAX 10000
#define ZEPS 1.0e-10
#define MOV3(a, b, c, d, e, f) (a)=(d);(b)=(e);(c)=(f);
double Simulator::dbrent(double ax, double bx, double cx, double tol, double* xmin)
{
    int iter, ok1, ok2;
    double a, b, d, d1, d2, du, dv, dw, dx, e = 0.0;
    double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = f1dim(x);
    dw = dv = dx = df1dim(x);

    for (iter = 1; iter <= ITMAX; iter++)
    {
        xm = 0.5 * (a + b);
        tol1 = tol * fabs(x) + ZEPS;
        tol2 = 2.0 * tol1;

        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
        {
            *xmin = x;
            return fx;
        }

        if (fabs(e) > tol1)
        {
            d1 = 2.0 * (b - a);
            d2 = d1;

            if (dw != dx)
            {
                d1 = (w - x) * dx / (dx - dw);
            }

            if (dv != dx)
            {
                d2 = (v - x) * dx / (dx - dv);
            }

            u1 = x + d1;
            u2 = x + d2;
            ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
            ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
            olde = e;
            e = d;

            if (ok1 || ok2)
            {
                if (ok1 && ok2)
                {
                    d = (fabs(d1) < fabs(d2) ? d1 : d2);
                }
                else if (ok1)
                {
                    d = d1;
                }
                else
                {
                    d = d2;
                }

                if (fabs(d) <= fabs(0.5 * olde))
                {
                    u = x + d;

                    if (u - a < tol2 || b - u < tol2)
                    {
                        d = SIGN2(tol1, xm - x);
                    }
                }
                else
                {
                    d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                }
            }
            else
            {
                d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
            }
        }
        else
        {
            d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
        }

        if (fabs(d) >= tol1)
        {
            u = x + d;
            fu = f1dim(u);
        }
        else
        {
            u = x + SIGN2(tol1, d);
            fu = f1dim(u);

            if (fu > fx)
            {
                *xmin = x;
                return fx;
            }
        }

        du = df1dim(u);

        if (fu <= fx)
        {
            if (u >= x)
            {
                a = x;
            }
            else
            {
                b = x;
            }

            MOV3(v, fv, dv, w, fw, dw)
            MOV3(w, fw, dw, x, fx, dx)
            MOV3(x, fx, dx, u, fu, du)
        }
        else
        {
            if (u < x)
            {
                a = u;
            }
            else
            {
                b = u;
            }

            if (fu <= fw || w == x)
            {
                MOV3(v, fv, dv, w, fw, dw)
                MOV3(w, fw, dw, u, fu, du)
            }
            else if (fu < fv || v == x || v == w)
            {
                MOV3(v, fv, dv, u, fu, du)
            }
        }
    }

    nrerror("Too many iterations in routine dbrent");
    return 0.0;
}
#undef ITMAX
#undef ZEPS
#undef MOV3

double Simulator::df1dim(double x)
{
    int j;
    double df1 = 0.0;
    double* xt, *df;

    xt = darray(ncom);
    df = darray(ncom);

    for (j = 0; j < ncom; j++)
    {
        xt[j] = pcom[j] + x * xicom[j];
    }

    dfunc(xt, df);

    for (j = 0; j < ncom; j++)
    {
        df1 += df[j] * xicom[j];
    }

    free_darray(df);
    free_darray(xt);
    return df1;
}

double Simulator::f1dim(double x)
{
    int j;
    double f, *xt;

    xt = darray(ncom);
    
    for (j = 0; j < ncom; j++)
    {
        xt[j] = pcom[j] + x * xicom[j];
    }
    
    f = func(xt);
    free_darray(xt);
    return f;
}

/*
 * ************************************
 * CONJUGATE GRADIENTS CODE ENDS HERE *
 * NOTHING SHALL GO BELOW THIS POINT  *
 * ************************************
 */