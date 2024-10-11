Collection of helper functions for simulations, mathematical operations using Eigen, building of block-structured matrices and much more...

## Logger
This is a `Logger` class used to collect a list of variables to log into csv files.
### Example
```cpp
#include <cpp-utils/simulation.h>

Eigen::MatrixXd m(2,3);
m << 1, 2, 3,
     4, 5, 6;
Eigen::Vector2d v = {2., 1.5};
double s = 3.14;
     
utils::Logger logger("data/");              // Insert local path where to create files
                                                  
logger.appendSuffix("toappend");            // String to append at the end: <FILENAME>_suffix.csv
logger.add(m, "mymatrix");
logger.add(v, "myvector", true);            // Transpose the vector (default is false)
logger.add(s, "myscalar", false, true);     // Add to the log list (default is true)

logger.logAll();

logger.logList();                           // Only log the variables added to the list

v(0) = 5;

logger.logAll();                            // All the variables are logged again, v(0) is changed


logger.log("myvector")                      // Log only one variable
```
