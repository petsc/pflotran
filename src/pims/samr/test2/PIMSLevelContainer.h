#include <vector>

using namespace std;

class PIMSLevelContainer{
public:
   PIMSLevelContainer(const int n);
   PIMSLevelContainer();
   ~PIMSLevelContainer();
   
   void *getPatchData(const int i){return patch_data[i];}

private:
   vector<void *> patch_data;

};
