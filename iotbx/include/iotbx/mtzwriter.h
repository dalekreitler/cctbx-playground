#ifndef MTZWRITER_H
#define MTZWRITER_H

#include <iostream>
#include <exception>
#include <string>

#include <sys/stat.h> // a cure for errors in irix compile
#include "cmtzlib.h"
#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/uctbx.h>
#include <iotbx/cppmtz.h>
#include <cctbx/sgtbx/space_group.h>

namespace af = scitbx::af;

namespace iotbx{ namespace mtz {

class MtzWriter {
private:
  CMtz::MTZ* mtz;
  CMtz::MTZXTAL* onextal;
  CMtz::MTZSET*  oneset;
public:
  MtzWriter();
  ~MtzWriter();

  void setTitle(const std::string&);
  void setSpaceGroup(const cctbx::sgtbx::space_group&);
  void oneCrystal(const std::string&,const std::string&,
                  const cctbx::uctbx::unit_cell&);
  void oneDataset(const std::string&,const double&);
  void write(const std::string&);
};
}} //namespaces

#endif /* MTZWRITER_H*/

