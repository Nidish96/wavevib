var colls = document.getElementsByClassName("collapsible");

function toggle_collapsible(coll){
    console.log(this);
    coll.classList.toggle("active");
	var content = coll.nextElementSibling;
    while(content){
	    if (content.style.display == "") {
	        content.style.display = "none";
	    } else {
	        content.style.display = "";
	    }
        content = content.nextElementSibling;
    }
}
for (let i = 0; i < colls.length; i++) {

    let coll = colls[i];
    toggle_collapsible(coll);
    coll.style.cursor = "pointer";
    coll.addEventListener("click", function() {
        toggle_collapsible(coll);

    });
}

var ols = document.getElementsByClassName("org-ol");

function toggle_ol_li_content(li){
    var children = li.children;
    if(children.length){
        for(j=1; j < children.length; j++){
            var content = children[j];
            if (content.style.display == "") {
	            content.style.display = "none";
	        } else {
	            content.style.display = "";
	        }
        }
    }
}

for (i = 0; i < ols.length; i++) {

    let ol = ols[i];
    var lis = ol.children;
    for(k=0;k<lis.length;k++){
        let li = lis[k];
        if(li.children.length > 1)
            li.style.cursor = "pointer";
        toggle_ol_li_content(li);
        li.addEventListener("click", function(e) {
            toggle_ol_li_content(li);
        });
    }

}
