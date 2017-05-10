#include "CellsNumber.h"

CellsNumber::CellsNumber(const char* name, const char* format) : Observer(name, format) {}

CellsNumber::CellsNumber(const CellsNumber& orig) : Observer(orig) {}

CellsNumber::~CellsNumber() {}

void CellsNumber::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double CellsNumber::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    return (double) cells.size();
}

DerivedRegister<CellsNumber> CellsNumber::reg("CellsNumber");

