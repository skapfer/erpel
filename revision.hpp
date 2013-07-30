#define REVISION_NUMBER_RAW "git-$Id$"
#include <string>

static std::string erpel_revision_number () {
    return "git-" + std::string (REVISION_NUMBER_RAW+9, 40);
}
