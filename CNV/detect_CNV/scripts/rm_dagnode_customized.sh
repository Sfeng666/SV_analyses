for id in $(cqd | awk '$2 ~ /1.1/{print $1}')
do
    echo $id
done

for id in $(cqd | awk '$6 == "R" {print}')
do
    echo $id
done