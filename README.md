# dtiQA_app
This includes everything required to build a docker and corresponding singularity container for the dtiQA ([topup_eddy_preprocess](https://github.com/justinblaber/topup_eddy_preprocess) + [dti_stats](https://github.com/justinblaber/dti_stats)) pipeline. 

# Build Instructions:
Just clone and run `build.sh`
```
git clone https://github.com/justinblaber/dtiQA_app.git
cd dtiQA_app/
./build.sh
```

# Run Instructions:
For docker:
```
sudo docker run --rm \
--runtime=nvidia \
-v $(pwd)/INPUTS/:/INPUTS/ \
-v $(pwd)/OUTPUTS:/OUTPUTS/ \
--user $(id -u):$(id -g) \
justinblaber/dtiqa
```
For singularity:
```
singularity run -e \
--nv \
-B INPUTS/:/INPUTS \
-B OUTPUTS/:/OUTPUTS \
dtiQA.img
```
For more detailed instructions on the inputs please visit: [https://justinblaber.org/dwmri-preprocessing/](https://justinblaber.org/dwmri-preprocessing/)
