#include "Simulator.h"

utils::Logger Simulator::simulator_logs("simulator");

Simulator::Simulator(const arguments& args) : number_of_cells(0), box(0, 0, 0),
    sb(args.render_file, args.surface_file, args.traj_file, args.stress_file),
    traj(args.traj_file, args.box_file), log_sim(args.output_file, args.ob_config_file)
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
    params.visc = args.visc;
    params.ttime = args.ttime;
    params.r_vertex = args.r_vertex;
    params.verlet_f = args.verlet_f;
    params.growth_rate = args.growth_rate;
    params.vc = args.vc;
    params.bud_d = args.bud_d;
    params.div_ratio = args.div_ratio;
    params.draw_box = args.draw_box;
    params.scale = args.scale_flag;
    params.dynamics = args.dynamics;
    params.nsteps = args.nsteps ? args.nsteps : (int)params.ttime / params.dt;
    params.platotype = args.platotype;
    params.v_disp_cut2 = params.r_vertex * params.r_vertex * (args.verlet_f - 1.0) * (args.verlet_f - 1.0);
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
    initCells(N, ra, rb, (char*)&"ms_kot");
}

void Simulator::initCells(int N, double ra, double rb, char* model_t)
{
    if (ra > rb)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Illegal arguments: ra > rb. Simulator will set: rb = ra \n";
        rb = ra;
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

//        if (Cell::membrane_test)
//        {
//            simulator_logs << utils::LogLevel::WARNING  << "MEMBRANE_TEST MODE\n";
//            double re = 1.0;
//            MembraneTriangulation mt;
//            tris = mt.triangulate(r0, re, 5);
//        }

        Cell newCell(tris);

        newCell.setEcc(params.E_cell);
        newCell.setNu(params.nu);
        newCell.setSpringConst(params.E_cell, params.th, params.nu, model_t);
        newCell.setBSprings(params.E_cell, params.th, params.nu);
        newCell.setDp(params.dp, params.ddp);
        newCell.setVertexR(params.r_vertex);
        newCell.setVerletR(params.verlet_f);
        newCell.setCellId(number_of_cells);
        newCell.setVisc(params.visc, params.dynamics);
        newCell.setInitR(r0);
        newCell.setGrowthRate(params.growth_rate);
        newCell.setBuddingVolume(params.vc);
        newCell.setBudDiameter(params.bud_d);
        newCell.setDivisionRatio(params.div_ratio);

//        if (Cell::membrane_test)
//        {
//            simulator_logs << utils::LogLevel::WARNING  << "MEMBRANE_TEST MODE\n";
//            double re = 1.0;
//            newCell._set_hooks(r0+re);
//        }

        pushCell(newCell);
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
    updateCells();

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

    bool resized = false;

    for (int i = 0; i < steps; i++)
    {

        if ( i % (steps / std::min(steps, 10) ) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * i / steps << "% OF THE SIMULATION IS DONE" "\n";
        }

        do
        {
            update_neighbors_list();
            integrate();
        }
        while ( check_min_force() );

        if ( (i + 1) % params.save_step == 0)
        {
            if (!params.scale)
            {
                traj.save(cells, getTotalVertices());
            }
            else
            {
                traj.save(cells, getTotalVertices(), box.getXmax() / box.getX(),
                          box.getYmax() / box.getY(), box.getZmax() / box.getZ());
            }

            traj.save_box(box, (i + 1) * params.dt);
        }

        if ( (i + 1) % params.log_step == 0)
        {
            log_sim.dumpState(box, cells);
        }

        resized = box.resize();

        if (resized)
        {
            domains.setBoxDim(box);
        }

        for (int i = 0; i < number_of_cells; i++)
        {
            cells[i].cellCycle(params.dt);
        }
    }

    log_sim.dumpState(box, cells);
    sb.saveStrainScript(cells, box);
    sb.saveStrainScript(cells, box);
    traj.close();
    log_sim.close();
}

void Simulator::calcForces()
{
//    double P0 = 0.025;
//    P0 = params.bud_d;
//    double R0 = 2.14;

//    double membrane_area = 3.14159265358979*R0*R0;
//    int num_in_R0 = cells[0].num_vertex(R0);
//    double force_per_vertex = P0 * membrane_area / num_in_R0;

//    double pulling_force = 0.01;
//    double hooks_pull_force = pulling_force / (cells[0].get_phooks_n() - 1);

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

//        if (Cell::membrane_test)
//        {
//            simulator_logs << utils::LogLevel::WARNING  << "MEMBRANE_TEST MODE\n";
//            for (int i = 0 ; i < number_of_cells; i++)
//            {
//                //cells[i].pull_vertex(pulling_force, 0.005);
//                //cells[i].pull_vertex(force_per_vertex, R0);
//                //cells[i].push_membrane(P0);
//                //cells[i].pull_membrane(hooks_pull_force);
//            }
//        }

        // CALCULATE INTER-CELLULAR FORCES
        #pragma omp for schedule(guided)

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
                //else if (params.nbhandler == 3)
                //{
                //    cells[i].calcNbForcesVL(cells[j], box);
                //}
                else
                {
                    cells[i].calcNbForcesON2(cells[j], box);
                }
            }
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

//        if (Cell::membrane_test)
//        {
//            simulator_logs << utils::LogLevel::WARNING  << "MEMBRANE_TEST MODE\n";
//            for (int i = 0 ; i < number_of_cells; i++)
//            {
//                //cells[i].voidForcesOutsideCircle(R0);
//                //cells[i].voidForcesForHooks();
//            }
//        }
    }
}

bool Simulator::verlet_condition()
{
    double disp = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            disp = cells[i].vertices[j].get_verlet_disp2();

            if (disp >= params.v_disp_cut2)
            {
                return true;
            }
        }

    }

    return false;
}

void Simulator::update_neighbors_list()
{
    if (params.nbhandler == 1)
    {
        if ( verlet_condition() )
        {
            rebuildVerletLists();
        }
    }
    else if (params.nbhandler == 2)
    {
        rebuildDomainsList();
    }

//    else if (params.nbhandler == 3)
//    {
//        if ( verlet_condition() )
//        {
//        }
//    }
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

//    else if (params.nbhandler == 3)
//    {
//    }

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

/*
 * INTEGRATORS
 */

void Simulator::integrate()
{
    (*this.*integrator)();
    updateCells();
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
            cells[i].vertices[j].r_c += dt * cells[i].vertices[j].f_c / visc;
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
            cells[i].vertices[j].r_p = cells[i].vertices[j].r_c;
            cells[i].vertices[j].f_p = cells[i].vertices[j].f_c;
        }
    }

    //move the whole time-step and calculate  forces
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].r_c += dt * cells[i].vertices[j].f_c / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].r_c = cells[i].vertices[j].r_p + 0.5 * dt * ( cells[i].vertices[j].f_p + cells[i].vertices[j].f_c) / visc;
        }
    }
}

void Simulator::midpointRungeKutta()
{
    double visc;
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
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].r_c += 0.5 * dt * cells[i].vertices[j].f_c / visc;
        }
    }

    calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            visc = cells[i].vertices[j].getVisc();
            cells[i].vertices[j].r_c = cells[i].vertices[j].r_p + dt * cells[i].vertices[j].f_c / visc;
        }
    }
}