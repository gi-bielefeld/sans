<!DOCTYPE html>

<head>
    <style>
        .headerBar 
        {
          position: fixed;
          width: 12%;
          left: 0%;
          top: 0%;
          z-index: 1;
          align: center;
          text-align: center;
          background-color: orange;
          display: block;
        }
        .tocHeader
        {
          margin: 0% auto 0% 17%;
          z-index: 1;
        }
        .tocExit
        {
          margin: 0% auto 0% 8.5%;
          z-index: 1;
          display: none;
    </style>
</head> 


<body>
<div id="headerBar" class="headerBar">
<div id="tocHeader" class="tocHeader" onclick="tocShow()">
Content &#11167;
</div>
<div id="tocExit" class="tocExit" onclick="tocHide()">
Content &#11165;
</div>
</div>

<script>
function showElement(id)
{
  document.getElementById(id).style.display = "block";
}
function hideElement(id)
{
  document.getElementById(id).style.display = "none";
}
function swapElement(oldId, newId)
{
  document.getElementById(oldId).style.display = "none";
  document.getElementById(newId).style.display = "block";
}
function tocShow()
{
  document.getElementById("tocHeader").style.display = "none";
  document.getElementById("tocExit").style.display = "block";
}

function tocHide()
{
  document.getElementById("tocExit").style.display = "none";
  document.getElementById("tocHeader").style.display = "block";
}

</script>

</body>