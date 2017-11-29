#include "Simulator.h"

utils::Logger Simulator::simulator_logs("simulator");
unsigned long Simulator::FORCE_EVALUATION_COUTER(0);

bool Simulator::RESTART_FLAG(false);

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
    catch (NotAllowedException& e)
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
    params.d = args.d;
    params.nbhandler = args.nb_flag;
    params.E_cell = args.E_cell;
    params.nu = args.nu;
    params.th = args.thickness;
    params.dt = args.dt;
    params.dp = args.dp;
    params.ddp = args.ddp;
    params.ttime = args.ttime;
    params.r_vertex = args.r_vertex;
    params.draw_box = args.draw_box;
    params.const_volume = args.const_volume;
    params.nsteps = args.nsteps ? args.nsteps : (int)params.ttime / params.dt;
    params.platotype = args.platotype;
    params.model_t = std::string(args.model_type);

    integrator = new Integrator(this, args.integrator_a);
    
    setTriangulator(args.tritype);
    box.setX(args.bsx);
    box.setY(args.bsy);
    box.setZ(args.bsz);
    box.setPbc(args.pbc);
    box.setEwall(args.E_wall);
    box.setNu(args.nu);
    box.setDefaultSchedule(params.nsteps, args.box_step, args.bsdx, args.bsdy, args.bsdz, 0.0, 0.0, 0.0);
    box.configureScheduler(args.sch_config_file);

    domains.setupDomainsList(getLengthScale( std::max(args.init_radius1, args.init_radius2) ), box);
    OsmoticForce::setVolumeFlag(args.osmotic_flag);
    OsmoticForce::setEpsilon(args.eps);
    Cell::no_bending = args.nobending;
    logParams();

}

Simulator::~Simulator()
{
    delete integrator;
}

void Simulator::diagnoseParams(arguments args)
{
    if (args.d < 0)
        throw NotImplementedException("NotImplementedException:\n"
                                      "Single point representation is not implemented yet. "
                                      "Simulator is about to terminate !");

    if ( (args.d == 0 && STRCMP(args.tritype, "rnd")) || (args.d == 0 && STRCMP(args.tritype, "plato")) )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Random and Plato triangulation cannot be used "
                                  "with single point representation.\n"
                                  "Simulation exits with EXIT_FAILURE status !");

    if ( (args.d == 0 && args.restart) || (args.d == 0 && args.analyze) )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Simulation for a single node representation"
                                  "cannot be restarted or post-analyzed!");

    if ( args.d == 0 && args.const_volume && args.nu != 0.5 )
        throw NotAllowedException("NotAllowedException:\n"
                                  "Single node representation:"
                                  " for const-volume option nu has to be equal to 0.5!");

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
    simulator_logs << utils::LogLevel::INFO  << "LOG_STEP="  << params.log_step << "\n";
    simulator_logs << utils::LogLevel::INFO  << "TRIANGULATOR="  << triangulator << "\n";
    simulator_logs << utils::LogLevel::FINE  << "TIME STEP(DT)="  << params.dt << " [s]\n";
    simulator_logs << utils::LogLevel::FINE  << "DEPTH="  << params.d << "\n";
    simulator_logs << utils::LogLevel::FINE  << "DP="  << params.dp << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "DDP="  << params.ddp << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "E CELL="  << params.E_cell << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "E BOX="  << box.getE() << " [MPa]\n";
    simulator_logs << utils::LogLevel::FINE  << "SURFACE_MODULUS="  << (params.E_cell * params.th) << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (CELL)="  << params.nu << "\n";
    simulator_logs << utils::LogLevel::FINE  << "POISSON'S_RATIO (BOX)="  << box.getNu() << "\n";
    simulator_logs << utils::LogLevel::FINE  << "R_VERTEX="  << params.r_vertex << " [micron]\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.PBC=" << (box.pbc ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINE  << "BOX.BOX_DRAW=" << (params.draw_box ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_FLAG=" << (OsmoticForce::getFlag() ? "true" : "false") << "\n";
    simulator_logs << utils::LogLevel::FINER << "OSMOTIC_EPS=" << OsmoticForce::getEpsilon() << "\n";
    simulator_logs << utils::LogLevel::FINER << "MAX_SCALE=" << domains.getMinLength() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.X="  << box.getX() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Y="  << box.getY() << "\n";
    simulator_logs << utils::LogLevel::FINER << "BOX.Z="  << box.getZ() << "\n";
}

void Simulator::initCells(int N, double r_min, double r_max, bool jam)
{
    if (r_min > r_max)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Illegal arguments: r_min > r_max. Simulator sets: r_min = r_max \n";
        r_max = r_min;
    }

    simulator_logs << utils::LogLevel::INFO  << "CELL MODEL: " << params.model_t << "\n";
    simulator_logs << utils::LogLevel::INFO  << "BENDING: " << (!Cell::no_bending ? "true" : "false") << "\n";

    double nx, ny, nz;
    bool flag = true;
    double rc = 2.0 * params.r_vertex;
    double rsum;
    double r0;

    while (number_of_cells < N)
    {
        r0 = uniform(r_min, r_max);
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
            addCell(r0);
            shiftCell(shift, number_of_cells - 1);
        }
    }

    if (jam)
    {
        simulator_logs << utils::LogLevel::INFO  << "SIMULATION STARTS FROM THE JAMMED PACKING\n";

        if (params.d == 0)
        {
            Packer::packCells(box, cells, params.th, false);
        }
        else
        {
            Packer::packCells(box, cells, params.th, true);
        }
    }

    if (params.d == 0)
    {
        simulator_logs << utils::LogLevel::WARNING  << "Single vertex representation\n";
        simulator_logs << utils::LogLevel::WARNING  << "Vertex radius reassignment: cell.vertex_r = cell.init_r\n";

        for (int i = 0; i < number_of_cells; i++)
        {
            cells[i].setVertexR( cells[i].getInitR() );
        }
    }

    restarter.saveTopologyFile(cells, params.model_t);
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

void Simulator::addCell(double r0)
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

        else if ( !triangulator.compare("rnd") )
        {
            RandomTriangulation rnd(10, 100, 0.1, 1000.0, params.r_vertex);
            tris = rnd.triangulate(r0);
        }

        Cell newCell(tris);
        newCell.setVertexR(params.r_vertex);
        newCell.setCellId(number_of_cells);
        newCell.setInitR(r0);

        newCell.setEcc(params.E_cell);
        newCell.setNu(params.nu);
        newCell.setSpringConst(params.E_cell, params.th, params.nu, params.model_t);
        newCell.setBSprings(params.E_cell, params.th, params.nu);
        newCell.setDp(params.dp, params.ddp);

        double radial_eps = 1.0 + (0.5 * (1 - params.nu)  * (newCell.getTurgor() * r0) / (params.E_cell * params.th));

        if (newCell.getNumberVertices() == 1)
        {
            radial_eps = 1.0;
        }

        newCell.setConstantVolume (radial_eps);

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
    simulator_logs << utils::LogLevel::INFO << "Simulation runs in [restart] mode. \n" ;
    restarter.registerVMap();
    restarter.readTopologyFile(cells);
    restarter.readLastFrame(cells);
    number_of_cells = cells.size();
    set_min_force();
    Simulator::RESTART_FLAG = true;
}

void Simulator::analyze()
{
    restarter.registerVMap();
    restarter.readTopologyFile(cells);
    number_of_cells = cells.size();
//
//    MIN_FORCE_SQ = 1e-8;
//    simulator_logs << utils::LogLevel::FINE  << "MIN_FORCE ARBITRARILY(in <<analyze>> mode) SET= "  << sqrt(MIN_FORCE_SQ) << " [units?]\n";

    uint frames_number = traj.countFramesNumber();
    simulator_logs << utils::LogLevel::INFO << " Number of frames in a trajectory file: " << (int) frames_number << "\n" ;

    std::vector<std::string> turgor_list = log_sim.readTurgorsFile();
    std::vector<std::string> boxsize_list = traj.read_saved_box();

    if (turgor_list.size() != frames_number)
    {
        simulator_logs << utils::LogLevel::SEVERE  << "CORRUPTED DATA: Turgor data doesn't match the number of frames. Observation analysis exists!";
        exit(EXIT_FAILURE);
    }

    if (boxsize_list.size() != frames_number)
    {
        simulator_logs << utils::LogLevel::SEVERE  << "CORRUPTED DATA: Box size data doesn't match the number of frames. Observation analysis exists!";
        exit(EXIT_FAILURE);
    }

    log_sim.registerObservers();
    log_sim.open();
    log_sim.printHeader();

    for (std::size_t i = 1; i <= frames_number; i++)
    {
        simulator_logs << utils::LogLevel::INFO << "[analyze] Processing frame number: " << i <<  "/" <<  (int) frames_number << "\n" ;
        restarter.readFrame(traj.getTrajFile(), cells, i);
        restarter.assignTurgors(turgor_list[i - 1], cells);
        restarter.assignBoxSize(boxsize_list[i - 1], box);
        updateCells();

        if ( i == 1 )
        {
            set_min_force();
            simulator_logs << utils::LogLevel::FINE  << "MIN_FORCE (in <<analyze>> mode) SET TO= "  << sqrt(MIN_FORCE_SQ) << " [units?]\n";
        }

        domains.setBoxDim(box);
        update_neighbors_list();
        log_sim.dumpState(box, cells, domains);
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
        rebuildDomainsList();
    }

    // LOGER READY TO WORK
    log_sim.registerObservers();
    log_sim.open();
    log_sim.printHeader();

    // TRAJECTORY FILE OPEND FOR DUMP
    traj.open();

    // IF SIMULATION RESTART - DON'T DUMP THE STATS
    if (!Simulator::RESTART_FLAG)
    {
        traj.save_traj(cells, getTotalVertices());
        log_sim.dumpState(box, cells, domains);
        saveTurgors();
        restarter.saveLastFrame(cells);
        restarter.saveTopologyFile(cells, params.model_t);
        traj.save_box(box, steps * params.dt);
        box.saveRemainingSchedule();
    }

    //================
    update_neighbors_list();
    calcForces();

    for (int i = 0; i < number_of_cells; i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            cells[i].vertices[j].f_p = cells[i].vertices[j].f_c;
        }
    }

    //===============

    bool resized = false;


    for (int step = 0; step < steps; step++)
    {
        if ( box.nthTodo() )
        {
            simulator_logs << utils::LogLevel::INFO << "The simulation has reached the end of the schedule @ the step " << step << "/" << steps << ".\n";
            break;
        }

        if ( step % (steps / std::min(steps, 10) ) == 0.0 )
        {
            simulator_logs << utils::LogLevel::INFO << 100.0 * step / steps << "% OF THE SIMULATION IS DONE" "\n";
        }

        unsigned long loop_couter = 0;

        do
        {
            do
            {
                integrate();
                loop_couter++;

                if ( (loop_couter + 1) % 5000 == 0)
                {
                    restarter.saveLastFrame(cells); // SO RESTARTING IS PRODUCTIVE
                }
            }
            while ( check_min_force() );
        }
        while ( check_const_volume() );
        
        integrator->resetParams(this);


        if ( (step + 1) % params.log_step == 0 )
        {
            // ** SAVE COORDINATES - i.e. "logging" coordinates
            traj.save_traj(cells, getTotalVertices());
            // **

            update_neighbors_list();
            log_sim.dumpState(box, cells, domains);
            saveTurgors();
            traj.save_box(box, (step + 1) * params.dt);
            restarter.saveLastFrame(cells);
            restarter.saveTopologyFile(cells, params.model_t);
        }

        if ( step < steps - 1 ) // DO NOT RESIZE ON THE LAST STEP
        {
            resized = box.resize( volumeFraction() );
        }
        else
        {
            resized = false;
        }

        if (resized)
        {
            domains.setBoxDim(box);
            box.saveRemainingSchedule();
        }

        recenterCells();

    }

    log_sim.dumpState(box, cells, domains);
    saveTurgors();
    traj.save_box(box, steps * params.dt);
    restarter.saveLastFrame(cells);
    traj.save_traj(cells, getTotalVertices());
    box.saveRemainingSchedule();

    traj.close();
    log_sim.close();

    simulator_logs << utils::LogLevel::FINEST << "Forces have been evaluated " << FORCE_EVALUATION_COUTER << " times.\n";
    simulator_logs << utils::LogLevel::FINEST << "Energy has been evaluated "  << Energy::ENERGY_EVALUATION_COUNTER << " times.\n";
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
        else if (params.nbhandler == 1)
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
    if (params.nbhandler == 1)
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

double Simulator::volumeFraction()
{
    double box_vol = box.getVolume();
    double cells_vol = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cells_vol += cells[i].calcVolume();
    }

    return (cells_vol / box_vol);

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

double Simulator::getLengthScale(double r_0)
{
    double maxscale = 2.0 * params.r_vertex;

    if (params.d == 0) // Vertex is a whole particle
    {
        maxscale = std::max(maxscale, 2.0 * r_0) ;
    }

    return maxscale;
}

/*
 * TRIANGULATION
 */
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
    else if (STRCMP (token, "rnd"))
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

    if (cells[0].getNumberVertices() == 1)
    {
        MIN_FORCE_SQ = 1e-12;
        MIN_FORCE_SQ = MIN_FORCE_SQ * MIN_FORCE_SQ;
    }

    Cell::FORCE_FRAC   = FORCE_FRAC;
    Cell::MIN_FORCE_SQ = MIN_FORCE_SQ;

    simulator_logs << utils::LogLevel::FINE  << "MIN_FORCE = "  << sqrt(MIN_FORCE_SQ) << " [units?]\n";
}

bool Simulator::check_min_force()
{
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
        double eps = 0.001; // 0.1% accuracy
        bool flag = false;

        for (int i = 0; i < number_of_cells; i++)
        {
            step = cells[i].checkVolumeCondition();

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

void Simulator::integrate()
{
    update_neighbors_list();
    integrator->integrate(this);
    updateCells();
}

void Simulator::saveTurgors()
{
    std::string turgorDumpFile = log_sim.getFileName() + ".turgor.out";

    double box_volume = box.getVolume();
    double box_x = box.getX();
    double box_y = box.getY();
    double box_z = box.getZ();
    double cells_vol = 0.0;

    FILE* os = fopen(turgorDumpFile.c_str(), "a");

    if ( os == NULL )
    {
        os = fopen(turgorDumpFile.c_str(), "w");
    }

    if ( os == NULL )
    {
        simulator_logs << utils::LogLevel::WARNING << "Can not open file:<<" << turgorDumpFile << "for writing.\n";
    }

    double turgor;
    Vector3D cm;

    for (uint i = 0; i < cells.size(); i++)
    {
        cells_vol += cells[i].calcVolume();
    }

    fprintf(os, "%-7.4f %5.3f %5.3f %5.3f %7.4f ", box_volume, box_x, box_y, box_z, cells_vol);


    for (uint i = 0; i < cells.size(); i++)
    {
        turgor = cells[i].getTurgor();
        cm = cells[i].getCm();

        fprintf(os, "%5.4f %6.4f %6.4f %6.4f ", turgor, cm.x, cm.y, cm.z);
    }

    fprintf(os, "%s", "\n");

    if ( os != NULL )
    {
        fclose(os);
    }
}

void Simulator::recenterCells()
{
    if (box.pbc)
    {
        for (uint i = 0; i < cells.size(); i++)
        {
            Vector3D shift = Box::recenteringVector( cells[i].getCm(), box );
            cells[i].addXYZ(shift);
        }
    }

    return;
}