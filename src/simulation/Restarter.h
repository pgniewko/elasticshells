#ifndef RESTARTER_H
#define	RESTARTER_H

#include <vector>
#include <string>

#include "Cell.h"
#include "utils/utils.h"

class Restarter {
public:
    Restarter(std::string, std::string);
    Restarter(const Restarter& orig);
    virtual ~Restarter();
    
    void saveTopologyFile(const std::vector<Cell>&) const;
    void readTopologyFile(const std::vector<Cell>&) const;
    
private:
    std::string topologyFile;
    std::string lastFrameFile;
    static utils::Logger restarter_logs;

};

#endif	/* RESTARTER_H */

