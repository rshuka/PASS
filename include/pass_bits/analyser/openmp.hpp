#pragma once

namespace pass
{
/**
 * Returns the evaluation time of a problem in nanoseconds
 *
 * NOTE: It is not a requirement to get the same output always because the CPU can
 * be less or more used by other processes running on the computer.
 * In the human mind, we can remember the solution of a math problem, though for a computer
 * the same process will always be something new, so, it is not required to get the same result always!
 *
 */
bool enable_openmp();

} // namespace pass
