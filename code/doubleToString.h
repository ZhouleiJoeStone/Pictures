#include <limits>
std::string DoubleToString(const double value, unsigned int precisionAfterPoint = 6)
{
   std::ostringstream out;
   // 清除默认精度
   out.precision(std::numeric_limits<double>::digits10);
   out << value;

   std::string res = std::move(out.str());
   auto pos = res.find('.');
   if (pos == std::string::npos)
       return res;

   auto splitLen = pos + 1 + precisionAfterPoint;
   if (res.size() <= splitLen)
       return res;

   return res.substr(0, splitLen);
}