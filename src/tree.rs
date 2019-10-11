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
    Node {
        data: T, 
        left_child: Box<Tree<T>>, 
        right_child: Box<Tree<T>>
    }
}

impl Tree<TreeNode> {
    fn from_linkage_matrix(l: Array2<f64>) -> Tree<TreeNode> {
        Tree::Empty
    }

    pub fn from_size(s: usize) -> Tree<TreeNode> {
        return if s > 0 {
            Tree::Node {
                data: TreeNode::from_usize(s),
                left_child: Box::new(Tree::Empty),
                right_child: Box::new(Tree::from_size(s - 1))
            }
        } else {
            Tree::Node {
                data: TreeNode::from_usize(s),
                left_child: Box::new(Tree::Empty),
                right_child: Box::new(Tree::Empty)
            }
        }

    }
    pub fn dfs(&self) -> Vec<TreeNode> {
        let mut result: Vec<TreeNode> = Vec::new();
        self._dfs(&mut result);
        return result
    }

    fn _dfs(&self, result: &mut Vec<TreeNode>) -> () {
        match self {
            Tree::Node { data, left_child, right_child } => {
                result.push(data.clone());
                left_child._dfs(result);
                right_child._dfs(result);
                println!("{:?}", data);
            },
            Tree::Empty => { /* Nothing */ }
        }
    }

    pub fn to_newick(&self) -> String {
        match self {
            Tree::Node { data, left_child, right_child } => {
                let is_tipnode = match **left_child {
                    Tree::Node {..} => false,
                    Tree::Empty => true
                };
                return if is_tipnode {
                    data.name.clone()
                } else {
                    format!(
                        "({},{})", 
                        left_child.to_newick(), 
                        right_child.to_newick()
                    ).to_string()
                }
            },
            Tree::Empty => {
                return "".to_string()
            }
        }
    }

}