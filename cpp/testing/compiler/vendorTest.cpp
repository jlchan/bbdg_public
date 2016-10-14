int main(int argc, char **argv){
#if defined(__clang__)
  return 1;

#elif defined(__ICC) || defined(__INTEL_COMPILER)
  return 2;

#elif defined(__GNUC__) || defined(__GNUG__)
  return 0;

#elif defined(__HP_cc) || defined(__HP_aCC)
  return 6;

#elif defined(__IBMC__) || defined(__IBMCPP__)
  return 4;

#elif defined(__PGI)
  return 5;

#elif defined(_CRAYC)
  return 8;

#elif defined(__PATHSCALE__) || defined(__PATHCC__)
  return 3;
#endif

  // Missing
  return 9;
}
