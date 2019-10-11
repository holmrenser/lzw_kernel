use ndarray::{Array2, Axis};

use crate::tree::{Tree, TreeNode};

pub trait NeighborJoining {
    fn from_neighbor_joining(
        mut dist: Array2<f64>,
        names: Vec<String>
    ) -> Tree<TreeNode>;
}

impl NeighborJoining for Tree<TreeNode> {
    fn from_neighbor_joining(
        mut dist: Array2<f64>,
        mut names: Vec<String>
    ) -> Tree<TreeNode> {
        let n = dist.shape()[0];
        let mut id = 0 as usize;
        let mut nodes: Vec<Tree<TreeNode>> = names.iter()
            .map(|n| {
                id += 1;
                Tree::Node(
                    TreeNode {
                        height: 1.0,
                        id: id,
                        name: n.to_string()
                    },
                    Box::new(Tree::Empty),
                    Box::new(Tree::Empty)
                )
            }).collect();
        loop {
            let q = make_q_mat(&dist);
            let pair = find_closest(&q);
            let (i1, i2) = pair;
            let node1 = nodes[i1].clone();
            let node2 = nodes[i2].clone();
            nodes.remove(i2);
            id += 1;
            nodes[i1] = Tree::Node(
                TreeNode {
                    height: 1.0,
                    id: id,
                    name: id.to_string()
                },
                Box::new(node1),
                Box::new(node2)
            );
            dist = make_new_dist(&dist, pair);
            if dist.len() <= 2 { break }
        }
        
        return nodes[0].clone()
    }
}


fn make_q_mat(dist: &Array2<f64>) -> Array2<f64> {
    let n = dist.len_of(Axis(0));
    let total_distance = dist.sum_axis(Axis(0));
    let _q: Vec<f64> = dist.indexed_iter()
        .map(|((row, col), d)| {
            return if row == col { 
                1.0 
            } else { 
                let scaled_dist = (n as f64 - 2.0) * d;
                (scaled_dist - &total_distance[row] - &total_distance[col])
            }
        }).collect();

    let q = Array2::from_shape_vec((n, n), _q).unwrap();
    return q
}

fn find_closest(q: &Array2<f64>) -> (usize, usize) {
    let ((row, col), _) = q.indexed_iter()
        .min_by(|(_, v1), (_, v2)| {
            v1.partial_cmp(v2).unwrap() // f64 forces the use of partial_cmp due to inf and nan
        }).unwrap();
    return if row < col { // return smallest index first, for use in make_new_dist
        (row, col) 
    } else {
        (col, row)
    }
}

fn make_new_dist(dist: &Array2<f64>, pair: (usize, usize)) -> Array2<f64> {
    let n = dist.len_of(Axis(0));
    let (node1, node2) = pair; // largest index last, use this to drop (f,g)
    let _d: Vec<f64> = (0..n.pow(2)).into_iter()
        .enumerate()
        .filter(|(i, _)| {
            let row = i / n;
            let col = i % n;
            return row != node2 && col != node2
        })
        .map(|(i, _)| {
            let row = i / n;
            let col = i % n;
            return if row == col {
                dist[[row, col]]
            } else if row == node1 {
                0.5 * (
                    dist[[node1, col]] 
                    + dist[[node2, col]] 
                    - dist[[node1, node2]]
                )
            } else if col == node1 {
                0.5 * (
                    dist[[node1, row]] 
                    + dist[[node2, row]] 
                    - dist[[node1, node2]]
                )
            } else {
                dist[[row, col]]
            }
        }).collect();

    let new_dist = Array2::from_shape_vec((n - 1, n - 1), _d).unwrap();
    
    return new_dist
}

pub fn neighbor_joining(mut dist: Array2<f64>) -> Array2<f64> {
    let n = dist.shape()[0];
    let mut node_names: Vec<usize> = (0..n).collect();
    let mut node_name = n;
    loop {
        let q = make_q_mat(&dist);
        let pair = find_closest(&q);
        let (node1, node2) = pair;
        let name1 = node_names[node1];
        let name2 = node_names[node2];
        node_names.remove(node2);
        node_names[node1] = node_name;
        node_name += 1;
        dist = make_new_dist(&dist, pair);
        if dist.len() <= 2 { break }
    }
    return dist
}