package uk.ac.shef.wit.geo.benchmark;

import org.openjdk.jmh.infra.BenchmarkParams;
import org.openjdk.jmh.infra.IterationParams;
import org.openjdk.jmh.profile.InternalProfiler;
import org.openjdk.jmh.results.AggregationPolicy;
import org.openjdk.jmh.results.IterationResult;
import org.openjdk.jmh.results.ScalarResult;

import java.util.ArrayList;
import java.util.Collection;

public class MaxMemoryProfiler implements InternalProfiler {

    @Override
    public String getDescription() {
        return "Max memory heap profiler";
    }

    @Override
    public void beforeIteration(BenchmarkParams benchmarkParams, IterationParams iterationParams) {

    }

    @Override
    public Collection<? extends ScalarResult> afterIteration(BenchmarkParams benchmarkParams,
                                                             IterationParams iterationParams,
                                                             IterationResult result) {

        Runtime runtime = Runtime.getRuntime();
        long memory = (runtime.totalMemory() - runtime.freeMemory())/1000000;

        Collection<ScalarResult> results = new ArrayList<>();
        results.add(new ScalarResult("memory min usage", memory, "MB", AggregationPolicy.MIN));
        results.add(new ScalarResult("memory max usage", memory, "MB", AggregationPolicy.MAX));
        results.add(new ScalarResult("memory avg usage", memory, "MB", AggregationPolicy.AVG));

        return results;
    }
}

