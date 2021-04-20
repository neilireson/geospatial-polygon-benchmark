## Benchmark of geospatial polygon indexes

[JMH](http://openjdk.java.net/projects/code-tools/jmh/) is used to compute the benchmarks

Currently, three libraries are included in the benchmark:

1. [Geotools](https://geotools.org)

2. [Lucene](https://lucene.apache.org)

2. [MongoDB](https://www.mongodb.com)

### Parameters

Include "-prof gc" for Memory performance

The benchmarks test both point and polygon intersection and vary

* Number of indexed polygons

### Results

On Geofabrik OSM England landuse data, querying with 10,000 random points/polygons.

Benchmark                                  Mode  Cnt   Score    Error  Units
GeotoolsBenchmark.pointIntersectsQuery    thrpt    3  13.640 ± 83.265  ops/s
GeotoolsBenchmark.polygonIntersectsQuery  thrpt    3   0.101 ±  0.422  ops/s
LuceneBenchmark.pointIntersectsQuery      thrpt    3   0.108 ±  0.514  ops/s
LuceneBenchmark.polygonIntersectsQuery    thrpt    3   0.092 ±  0.117  ops/s
MongoDbBenchmark.pointQuery               thrpt    3   0.095 ±  0.049  ops/s
MongoDbBenchmark.polygonQuery             thrpt    3   0.028 ±  0.022  ops/s
PostgisBenchmark.pointQuery               thrpt    3   0.091 ±  0.142  ops/s
PostgisBenchmark.polygonQuery             thrpt    3   0.065 ±  0.068  ops/s