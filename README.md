# DRTS
## Configuration

### Required Configuration

You **must configure** the following variable in the `int main()` function before running the program:

- `string filename`: Specifies the benchmark instance to be tested.

### Optional Configuration

The following variables can be configured optionally. If not specified, their default values will be used.

- `double time_limit_seconds`: Maximum runtime limit in seconds. Default: `1000`.
- `int fixed_maxIter`: Maximum number of iterations. Default: `2000000`.
- `int th`: Number of consecutive iterations without improving the best-known solution before triggering perturbation. Default: `200`.
- `double min_k1`: Range of the proportion of vertices removed from the dominating set during perturbation. Default range: `0.05–0.15`.
- `double min_k2`: Range of the proportion of vertices added to the dominating set from outside the current dominating set during perturbation. Default range: `0.30–0.50`.
- `global_seed`: Random seed. By default, a random seed is generated.

## Output Explanation

During execution, the program outputs information about the current iteration, the operation performed, the resulting minimum spanning tree (MST), and the best solution found so far.

**Example output:**

```text
Iter 5218 add 82 -> MST 851.69  best_result 751.75  time to best 23.2652  size 8  best solution number : 2
```

The fields are interpreted as follows:

- `Iter 5218`: The current iteration number is 5218.
- `add 82`: Vertex 82 is added to the dominating set.
- `MST 851.69`: The weight of the minimum spanning tree induced by the updated dominating set is 851.69.
- `best_result 751.75`: The best objective value found so far is 751.75.
- `time to best 23.2652`: The best objective value of 751.75 was first found after 23.2652 seconds of execution.
- `size 8`: The corresponding best solution contains 8 vertices in its dominating set.
- `best solution number : 2`: Two distinct solutions with an objective value of 751.75 have been found.

