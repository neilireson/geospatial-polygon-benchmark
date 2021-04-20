## Benchmark of geospatial polygon indexes

[JMH](http://openjdk.java.net/projects/code-tools/jmh/) is used to compute the benchmarks

Currently, three libraries are included in the benchmark:

1. [Geotools](https://geotools.org)

2. [Lucene](https://lucene.apache.org)

2. [MongoDB](https://www.mongodb.com)

### Parameters

The benchmarks test both point and polygon intersection and vary

* Number of indexed polygons

### Results

Benchmark                                      Mode      Score   Error   Units
GeotoolsBenchmark.pointIntersectsQuery        thrpt     12.849           ops/s
GeotoolsBenchmark.polygonIntersectsQuery      thrpt      0.073           ops/s
LuceneBenchmark.pointIntersectsQuery          thrpt      0.073           ops/s
LuceneBenchmark.polygonIntersectsQuery        thrpt      0.085           ops/s
MongoDbBenchmark.pointQuery                   thrpt      0.088           ops/s
MongoDbBenchmark.polygonQuery                 thrpt      0.030           ops/s
PostgisBenchmark.pointQuery                   thrpt      0.090           ops/s
PostgisBenchmark.polygonQuery                 thrpt      0.072           ops/s
