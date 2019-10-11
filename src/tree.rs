use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct TreeNode {
    pub height: f64,
    pub id: usize,
    pub name: String
}

impl TreeNode {
    fn from_usize(s: usize) -> TreeNode {
        return TreeNode {
            height: s as f64,
            id: s,
            name: "".to_string()
        }
    }
}

#[derive(Debug, Clone)]
pub enum Tree<T> {
    Empty,
    Node(T, Box<Tree<T>>, Box<Tree<T>>)
}

impl Tree<TreeNode> {
    fn from_linkage_matrix(l: Array2<f64>) -> Tree<TreeNode> {
        Tree::Empty
    }

    pub fn from_size(s: usize) -> Tree<TreeNode> {
        return if s > 0 {
            Tree::Node(
                TreeNode::from_usize(s),
                Box::new(Tree::Empty),
                Box::new(Tree::from_size(s - 1))
            )
        } else {
            Tree::Node(
                TreeNode::from_usize(s),
                Box::new(Tree::Empty),
                Box::new(Tree::Empty)
            )
        }

    }
    pub fn dfs(&self) -> Vec<TreeNode> {
        let mut result: Vec<TreeNode> = Vec::new();
        self._dfs(&mut result);
        return result
    }

    fn _dfs(&self, result: &mut Vec<TreeNode>) -> () {
        match self {
            Tree::Node(val, left_child, right_child) => {
                result.push(val.clone());
                left_child._dfs(result);
                right_child._dfs(result);
                println!("{:?}", val);
            },
            Tree::Empty => { /* Nothing */ }
        }
    }

    fn to_newick(&self) -> String {
        return "(1,2);".to_string();
    }

}