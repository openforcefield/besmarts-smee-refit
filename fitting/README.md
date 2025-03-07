# Fitting

The notebook here `example-smee-descent.ipynb` shows quickly how to load a dataset and what an entry might look like.

For how to fit a force field to optimization data from a SMIRNOFF force field, here's an [example I put together for the IRL Irvine meeting](https://openforcefield.atlassian.net/wiki/spaces/MEET/pages/3440508935/Hackathon+How+to+train+your+force+field+with+smee) (`run-smee-fit-from-qca-data-commented.ipynb` where you can largely follow on from the "Assign parameters to molecules in the dataset" heading.

The only note is that the `descent.targets.energy.predict` function would have to be re-written to not include forces in the objective and prediction if they're not in the data.