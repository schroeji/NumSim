#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>


#include <assert.h>
#include "../src/typedef.hpp"
#include <tuple>
#include <vector>





std::vector< real_t >
convertToReal_t
(
   std::string& stringToConvert   
)
{
    std::vector< real_t  > output;
    //std::index_t pos = stringToConvert.find( ',' );
    std::string temp;
    for( auto it : stringToConvert )
    {
       if( it == ','  )
       {
          output.push_back( std::stod( temp ) );
          temp.clear();
       }
       else
       {
          temp.push_back( it );   
       }
    }
    output.push_back( std::stod( temp ) );
    return output;
}



std::vector< std::tuple< real_t, real_t, real_t > >
load
(
   const std::string path,
   const real_t      re
)
{
  std::fstream file( path , std::ios::in /*| std::ios::nocreate*/ );
  
  assert( file.good() );
  std::vector< std::tuple< real_t, real_t, real_t > > output;
  
  std::string line;
  std::vector< real_t > temp;
  while ( std::getline( file, line ) ) 
  {
      temp = convertToReal_t( line );
      if( 1 == temp.size() )
      {
          assert( temp.at(0) - re < 1.0e-10 );
      }
      else if( 3 == temp.size() )
      {
         output.push_back( std::make_tuple( temp.at( 0), temp.at( 1 ), temp.at( 2 ) ) );
      }
      else
      {
         assert( false );   
      }
  }
  file.close( );
  return output;

}



void 
computeMean
(
   const real_t h    
)
{
    std::vector< std::tuple< real_t, real_t, real_t > >  temp;
    real_t erw = 0.0;
    std::vector< real_t > firstCellMean;
    std::vector< real_t > secondCellMean;
    std::vector< real_t > thirdCellMean;
    for( real_t re = 1000.000000; re <= 2000.0 ; re +=h )
    {
       std::string path = "uvalues_uniform/re_" + std::to_string( re );
       temp = load( path, re );
       real_t weight = h;
       if( re == 1000.0 || re == 2000.0 )
       {
           weight *= 0.5;
           //std::cout << weight << std::endl;
           if( firstCellMean.size() != temp.size() )
           {
               for( index_t i = 0; i < temp.size(); i++ )
               {
                  firstCellMean.push_back( 0.0 );
                  secondCellMean.push_back( 0.0 );
                  thirdCellMean.push_back( 0.0 );
               }
           }
       }
       
       for( index_t i = 0; i < temp.size(); i++ )
       {
          firstCellMean.at( i ) += weight * std::get< 0 >( temp.at( i ) );
          secondCellMean.at( i ) += weight * std::get< 1 >( temp.at( i ) );
          thirdCellMean.at( i ) += weight * std::get< 2 >( temp.at( i ) );
       }
    }
    
    std::fstream f;
    f.open("Mean" + std::to_string( h ), std::ios::out | std::ios::app);
    for( index_t i = 0; i < firstCellMean.size(); i++ )
    {
    f << firstCellMean.at( i ) * 1.0e-3 << "," << secondCellMean.at( i )* 1.0e-3 << "," << thirdCellMean.at( i )* 1.0e-3 << std::endl;
    }
    f.close();
    
}




void
computeStandardDeviation
(
   const real_t h
)
{
    std::vector< std::tuple< real_t, real_t, real_t > >  temp;
    std::vector< std::tuple< real_t, real_t, real_t > > mean = load( "Mean" + std::to_string(h), 0.0 );
    std::vector< real_t > firstCellStandardDeviation;
    std::vector< real_t > secondCellStandardDeviation;
    std::vector< real_t > thirdCellStandardDeviation;
    for( real_t re = 1000.000000; re <= 2000.0 ; re +=h )
    {
       std::string path = "uvalues_uniform/re_" + std::to_string( re );
       //std::cout << path << std::endl;
       temp = load( path, re );
       real_t weight = h;
       if( re == 1000.0 || re == 2000.0 )
       {
           weight *= 0.5;
           if( firstCellStandardDeviation.size() != temp.size() )
           {
               for( index_t i = 0; i < temp.size(); i++ )
               {
                  firstCellStandardDeviation.push_back( 0.0 );//resize( temp.size() );
                  secondCellStandardDeviation.push_back( 0.0 );//resize( temp.size() );
                  thirdCellStandardDeviation.push_back( 0.0 );//resize( temp.size() );   
               }
           }
       }
       for( index_t i = 0; i < temp.size(); i++ )
       {
          firstCellStandardDeviation.at( i ) += weight 
                                                * ( std::get< 0 >( temp.at( i ) ) - std::get< 0 >( mean.at( i ) ) ) 
                                                * ( std::get< 0 >( temp.at( i ) ) - std::get< 0 >( mean.at( i ) ) );
          secondCellStandardDeviation.at( i ) += weight 
                                                * ( std::get< 1 >( temp.at( i ) ) - std::get< 1 >( mean.at(i) ) )
                                                * ( std::get< 1 >( temp.at( i ) ) - std::get< 1 >( mean.at(i) ) );
          thirdCellStandardDeviation.at( i ) += weight 
                                                * ( std::get< 2 >( temp.at( i ) ) - std::get< 2 >( mean.at( i ) ) )
                                                * ( std::get< 2 >( temp.at( i ) ) - std::get< 2 >( mean.at( i ) ) );
       }
    }
    
    std::ofstream f;
    f.open("StandardDeviation" + std::to_string( h ), std::ios::out | std::ios::app);
    for( index_t i = 0; i < firstCellStandardDeviation.size(); i++ )
    {
        f << firstCellStandardDeviation.at( i )* 1.0e-3 
          << "," 
          << secondCellStandardDeviation.at( i )* 1.0e-3 
          << "," 
          << thirdCellStandardDeviation.at( i ) * 1.0e-3
          << std::endl;
    }
    f.close();
}




int main( int argc, char **argv )
{
    for( real_t h : {20.0, 10.0, 5.0 } )
    {
       computeMean( h );
       computeStandardDeviation( h );
    }

   return 0;
}



