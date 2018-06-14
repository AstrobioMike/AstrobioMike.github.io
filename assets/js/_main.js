// Scrolling changes background-color of li element in table-of-contents
$(document).ready(function() {
  var screenHeight = window.innerHeight || document.documentElement.clientHeight || document.body.clientHeight;

  window.addEventListener('scroll', function(event) {
    if (!$('#table-of-contents').length)
      return;

    var min_offset = null;
    var target_id = null;

    $(':header').each(function (index, item) {
      if ($(item).attr('id')) {
        var offset = $(item).offset().top - $(document).scrollTop(); 

        if (offset < (screenHeight / 2)) {
          if (min_offset == null || Math.abs(offset) < Math.abs(min_offset)) {
            min_offset = offset;
            target_id = $(item).attr('id');
          }
        }
      };
    });
    $("#markdown-toc a").removeClass('header-active');

    if (target_id) {
      $("#markdown-toc a[href$='#"+target_id+"']").addClass('header-active');
    }
  }, true);

  // this is for smooth scrolling
  $('#markdown-toc a[href^="#"]').on('click', function(event) {
      var anchor = this;
      var target = $(anchor.getAttribute('href'));
      if( target.length ) {
          event.preventDefault();
          $('html, body').stop().animate({
              scrollTop: target.offset().top
          }, 500)
          .promise()
          .done(function(){ location.href = anchor.getAttribute('href'); });
      }
  });
});
