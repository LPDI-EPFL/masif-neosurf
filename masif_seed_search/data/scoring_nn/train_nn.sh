source /work/upcorreia/bin/load_masif_environment_gpu.sh
masif_root=../../../
masif_seed_root=../../
masif_source=$masif_root/masif/source/
masif_seed_source=$masif_seed_root/source/
masif_data=$masif_root/masif/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_seed_source
python3 $masif_seed_source/train_alignment_evaluation_nn.py 12A 0_1_2_3
#python3 $masif_seed_source/train_alignment_evaluation_nn.py 12A $1
